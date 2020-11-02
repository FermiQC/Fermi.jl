using Combinatorics
using SparseArrays
using TensorOperations
using LinearAlgebra
using ArnoldiMethod

function CASCI{T}(Alg::xACI) where T <: AbstractFloat
    # we need this implementation
    @output "Getting molecule...\n"
    #molecule = Molecule()
    @output "Computing AO Integrals...\n"
    #aoint = ConventionalAOIntegrals()

    @output "Calling RHF module...\n"
    refwfn = Fermi.HartreeFock.RHF()
    CASCI{T}(refwfn, Alg)
end

function CASCI{T}(ci::CASCI, Alg::xACI) where T <: AbstractFloat
    @output "Using previous CASCI wave function as starting point\n"
    CASCI{T}(ci.ref, Alg, ci=ci)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, ci::CASCI, Alg::xACI) where T <: AbstractFloat
    @output "Using previous CASCI wave function as starting point\n"
    CASCI{T}(refwfn, Alg, ci=ci)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, Alg::xACI; ci = nothing) where T <: AbstractFloat
    @output "Generating Integrals for CAS computation...\n"
    #aoint = ConventionalAOIntegrals()
    ints = refwfn.ints
    @output "Transforming Integrals for CAS computation...\n"
    frozen = Fermi.CurrentOptions["cas_frozen"]

    nmo = refwfn.ndocc + refwfn.nvir

    act_elec = 2*(refwfn.ndocc - frozen)

    if act_elec < 0
        error("\nInvalid number of frozen orbitals ($frozen) for $(2*refwfn.ndocc) electrons.")
    end

    # Active = -1 means FCI, with frozen
    if Fermi.CurrentOptions["cas_active"] == -1
        active = nmo - frozen
    else
        active = Fermi.CurrentOptions["cas_active"]
    end

    if active ≤ act_elec/2
        error("\nNumber of active orbitals ($active) too small for $(act_elec) active electrons")
    end

    if active+frozen > nmo
        error("\nNumber of active ($active) and frozen orbitals ($frozen) greater than number of orbitals ($nmo)")
    end

    h = Fermi.Integrals.transform_fock(ints["T"] + ints["V"], ints.orbs["FU"], ints.orbs["FU"])
    V = Fermi.Integrals.transform_eri(ints["μ"], ints.orbs["FU"], ints.orbs["FU"], ints.orbs["FU"], ints.orbs["FU"])

    aoint = nothing

    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg, ci=ci)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::xACI; ci=nothing) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    ttotal = @elapsed begin
    @output "\n    • Computing CI with the extended ACI algorithm.\n\n"
    act_range = (frozen+1):(active+frozen)
    σ = Fermi.CurrentOptions["σ"]
    γ = Fermi.CurrentOptions["γ"]
    ζ = Fermi.CurrentOptions["ζ"]
    ζsize = Fermi.CurrentOptions["ζsize"]
    pe = Fermi.CurrentOptions["aci_print_screen"]
    nmo = refwfn.ndocc + refwfn.nvir

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:           {:3d}\n" frozen
    @output "Active Electrons:          {:3d}\n" act_elec
    @output "Active Orbitals:           {:3d}\n" active
    @output "Total number of Orbitals:  {:3d}\n" nmo

    # Start reference space as HF
    zeroth = repeat('1', frozen)*repeat('1', Int(act_elec/2))
    if ci == nothing
        P = [Determinant(zeroth, zeroth)]
        Pcoef = [1.0]
    else
        P, Pcoef = coarse_grain(ci.dets, ci.coef, γ, σ)
    end
    E = refwfn.energy - refwfn.molecule.Vnuc
    ΔE = 1.0
    ite = 1

    @output repeat("=",50)*"\n"
    Nα = Int(act_elec/2)
    Nβ = Int(act_elec/2)
    Lenny = length(P)
    M = nothing
    ϵsum = nothing
    ϵest = nothing
    oldP = nothing
    cflag = false
    while true
        @output " → Iteration {}\n\n" ite
        @output "   • P\n"
        @output "Initial model space size: {}\n\n" length(P)

        @output "   • P ⇒ F\n"
        @output "Generating First Order Interacting Space...\n"
        t = @elapsed F = get_fois(P, Nα, Nβ, act_range)
        @output "FOIS size:                {}\n" length(F)
        @output "FOIS contructed in {:5.5f} s.\n\n" t

        @output "   • F ⇒ Q\n"
        @output "Screening FOIS using 2-D Hamiltonian\n" σ
        t = @elapsed Fe = ϵI(F, P, Pcoef, E, h, V, act_elec, active)
        @output "Screen complete in {:5.5} s.\n" t
        @output "Sorting F space...\n"
        _t = @elapsed begin
            Fperm = zeros(Int, length(Fe))
            sortperm!(Fperm, Fe, by=abs)
            reverse!(Fperm)
            Fe = Fe[Fperm]
            F = F[Fperm]
        end
        @output "Sorted in {:5.5f} s.\n" _t
        @output "Filtering F..."
        t = @elapsed begin
        ϵest = 0.0
        ϵsum = 0.0
        while true
            if length(Fe) == 0 
                #then no determinants were deemed important - exit ACI
                break
            end
            if ϵsum ≤ σ
                ϵest += Fe[end]
                ϵsum += abs(Fe[end])
                pop!(Fe)
                pop!(F)
            else
                ϵest -= Fe[end]
                ϵsum -= abs(Fe[end])
                break
            end
        end
        end
        @output " Secondary space (Q) built in {:5.5f}\n\n" t
        @output "Size of Q {}\n" length(Fe)

        @output "   • M = P ∪ Q\n"
        ΔE = -E
        M = vcat(P, F)
        @output "Model space size: {}\n" length(M)
        @output "Updating model space...\n" length(M)
        t = @elapsed E, Pcoef, P = update_model_space(M, h, V)
        @output " Model space updated in {:5.5f}\n" t
        ΔE += E
        @output "Model Space Energy           {:15.10f}\n" E + refwfn.molecule.Vnuc
        @output "Energy Change                {:15.10f}\n" ΔE

        if oldP == Set(P) 
            cflag = true
            break
        end
        ite += 1
        if ite > 30
            break
        end
        oldP = Set(deepcopy(P))
        P, Pcoef = coarse_grain(P, Pcoef, γ, σ)

        @output "Final coarse grained model space size is {}\n" length(P)
        @output repeat("=",50)*"\n"
    end
    end #@elapsed

    @output repeat("=",50)*"\n"
    if cflag
        @output "🔥🔥🔥🔥🔥 ACI procedure has converged. 🔥🔥🔥🔥🔥\n"
    else
        @output "😲😲😲😲😲 ACI procedure has failed!!!! 😲😲😲😲😲\n"
    end
    @output "Computation finished in {:5.5} seconds.\n" ttotal
    @output "Model space size: {}\n" length(M)
    @output "E[ACI:{}]     = {:15.10f}\n" σ E + refwfn.molecule.Vnuc
    @output "E[ACI:{}]+PT2 = {:15.10f}\n" σ E + refwfn.molecule.Vnuc + ϵest
    @output repeat("=",51)*"\n\n"

    
    ζdets = Determinant[]
    if ζsize == nothing
        @output "\nSelecting reference determinants using ζ = {:3.2f}" ζ
    else
        @output "\nReferemce determinants are going to be the {} most important ones" ζsize
    end

    acum = 0.0
    for i in eachindex(P)
        push!(ζdets, P[i])
        acum += Pcoef[i]^2
        if acum > ζ && (ζsize==nothing) > ζ
            break
        end

        if length(ζdets) == ζsize
            break
        end
    end

    # Get top determinants from model space (based on parameter ζ)
    @output "\n • {:d} ζ-determinants:\n\n" length(ζdets)

    for i in 1:(length(ζdets))
        @output "{:15.5f}      {}\n" Pcoef[i]  detstring(P[i], frozen+active)
    end
    @output "\n∑ C^2 = {:3.2f}\n" acum

    # Generate all single and double excitations from these selected determinants. Excitations are taken in the full orbital space.
    full_range = (frozen+1):size(h,1)
    println(full_range)
    ζSD = get_fois(ζdets, Nα, Nβ, (frozen+1):nmo)

    # Remove determinants that are already in the model space
    filter!(d->!(d in P), ζSD) 

    @output "Size of the ζ-FOIS space: {:5.5d}\n" length(ζSD)

    for d in ζSD
        if sum(αlist(d)) != Nα
            println("Invalid det found")
            showdet(d)
            @assert 1==2
        elseif sum(βlist(d)) != Nβ
            println("Invalid det found")
            showdet(d)
            @assert 1==2
        end
    end

    @output "Screening ζ-FOIS using 2-D Hamiltonian\n" σ
    t = @elapsed Fe = ϵI(ζSD, P, Pcoef, E, h, V)
    @output "Screen complete in {:5.5} s.\n" t
    @output "Sorting ζ-determinants space...\n"
    _t = @elapsed begin
        Fperm = zeros(Int, length(Fe))
        sortperm!(Fperm, Fe, by=abs)
        reverse!(Fperm)
        Fe = Fe[Fperm]
        ζSD = ζSD[Fperm]
    end
    @output "Sorted in {:5.5f} s.\n" _t
    @output "Filtering ζ-determinants...\n"
    t = @elapsed begin
    ϵest = 0.0
    ϵsum = 0.0
    while true
        if length(ζSD) == 0 
            #then no determinants were deemed important - exit ACI
            break
        end
        if ϵsum ≤ σ
            ϵest += Fe[end]
            ϵsum += abs(Fe[end])
            pop!(Fe)
            pop!(ζSD)
        else
            ϵest -= Fe[end]
            ϵsum -= abs(Fe[end])
            break
        end
    end
    end
    @output " Expanded space built in {:5.5f}\n\n" t
    @output "Size of filtered ζ-FOIS {}\n" length(Fe)

    M = vcat(P, ζSD)
    @output "Final extended space size: {}\n" length(M)
    @output "Performing final diagonalization...\n"
    t = @elapsed E, Pcoef, P = update_model_space(M, h, V, complete=false)
    @output "Final xACI({}) Energy           {:15.10f}\n" σ E + refwfn.molecule.Vnuc
    @output "Final xACI({}) Energy+PT2       {:15.10f}\n" σ E + refwfn.molecule.Vnuc + ϵest
    
    CASCI{T}(refwfn, E+refwfn.molecule.Vnuc, P, Pcoef)
end

@fastmath @inbounds function get_fois(dets::Array{Determinant,1}, Nα::Int, Nβ::Int, act_range::UnitRange{Int64})::Array{Determinant,1}

    # Check for the needed type of Int
    if length(act_range) < 60
        one = Int64(1)
    else
        one = Int128(1)
    end

    # Ns must be > 1
    αoccs = [zeros(Int,Nα) for i=1:Threads.nthreads()]
    βoccs = [zeros(Int,Nβ) for i=1:Threads.nthreads()]
    αunos = [zeros(Int,length(act_range)-Nα) for i=1:Threads.nthreads()]
    βunos = [zeros(Int,length(act_range)-Nβ) for i=1:Threads.nthreads()]

    lf_per_det = (length(αoccs[1])^2*length(αunos[1])^2 + length(αoccs[1])*length(αunos[1])
                       + length(βoccs[1])^2*length(βunos[1])^2 + length(βoccs[1])*length(βunos[1])
                       + length(αoccs[1])*length(αunos[1])*length(βoccs[1])*length(βunos[1]))
    lf_crit = Int(round(length(dets)*lf_per_det))
    fois = [Determinant(0,0) for i=1:lf_crit]
    @sync for _DI in eachindex(dets)
        Threads.@spawn begin
            d = dets[_DI]
            DI = (_DI-1)*lf_per_det + 1
            ct = 0
            id = Threads.threadid()
            αocc = αoccs[id]
            βocc = βoccs[id]
            αuno = αunos[id]
            βuno = βunos[id]
            αocc!(d, act_range, αocc)
            βocc!(d, act_range, βocc)
            αvir!(d, act_range, αuno)
            βvir!(d, act_range, βuno)
            # Get αα -> αα excitations
            for i in αocc
                for a in αuno
                    newα = (d.α ⊻ (one<<(i-1))) | (one<<(a-1)) 
                    _det = Determinant(newα, d.β)
                    fois[DI+ct] = _det
                    ct += 1
                    for j in αocc
                        if j ≥ i
                            break
                        end
                        for b in αuno
                            if b ≥ a
                                break
                            end
                            newestα = (newα ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(newestα, d.β)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
            # Get ββ -> ββ excitations
            for i in βocc
                for a in βuno
                    newβ = (d.β ⊻ (one<<(i-1))) | (one<<(a-1)) 
                    _det = Determinant(d.α, newβ)
                    fois[DI+ct] = _det
                    ct += 1
                    for j in βocc
                        if j ≥ i
                            break
                        end
                        for b in βuno
                            if b ≥ a
                                break
                            end
                            newestβ = (newβ ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(d.α, newestβ)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
            # Get αβ -> αβ excitations
            for i in αocc
                for a in αuno
                    for j in βocc
                        for b in βuno
                            newα = (d.α ⊻ (one<<(i-1))) | (one<<(a-1)) 
                            newβ = (d.β ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(newα, newβ)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
        end #Threads.@spawn 
    end
    fois = filter((x)->x != Determinant(0,0), fois)
    fois = Set(fois)
    setdiff!(fois, dets)
    fois = collect(fois)
    return fois end

#lexicographic premutations generation, By Donald Knuth
function lpermutations(a::BitArray)
  b=BitArray[]
  sort!(a)
  n=length(a)
  while(true)
    push!(b,copy(a))
    j=n-1
    while(a[j]>=a[j+1])
      j-=1
      j==0 && return(b)
    end
    l=n
    while(a[j]>=a[l])
      l-=1
    end
    tmp=a[l]
    a[l]=a[j]
    a[j]=tmp
    k=j+1
    l=n
    while(k<l)
      tmp=a[k]
      a[k]=a[l]
      a[l]=tmp
      k+=1
      l-=1
    end
  end
end

function ϵI(Fdets, P::Array{Determinant,1}, Pcoef::Array{Float64,1}, Ep::T, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat
    @output "Starting ϵI...\n"
    Fe = zeros(length(Fdets))
    N = sum(αlist(P[1]))
    αinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    βinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    #αind = Array{Int64,1}(undef,N) 
    #βind = Array{Int64,1}(undef,N) 
    @sync for i in eachindex(Fdets)
    #for i in eachindex(Fdets)
        begin
        D1 = Fdets[i]
        id = Threads.threadid()
        αind = αinds[id]
        βind = βinds[id]
        αindex!(D1, αind)
        βindex!(D1, βind)
        Ei = Hd0(αind, βind, h, V)
        Δ = Ei - Ep
        Vint = 0.0
        for j in eachindex(P)
            D2 = P[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                continue 
            elseif el == 2
                Vint += Pcoef[j]*Hd2(D1, D2, V, αexc)
            elseif el == 1
                Vint += Pcoef[j]*Hd1(αind, βind, D1, D2, h, V, αexc)
            end
        end
        
        @fastmath Fe[i] = Δ/2 - √((Δ^2)/4 + Vint^2)
        end #Threads.@spawn
    end
    return Fe
end
