using Combinatorics
using SparseArrays
using TensorOperations
using LinearAlgebra
using ArnoldiMethod

function CASCI{T}(Alg::ACI) where T <: AbstractFloat
    # we need this implementation
    @output "Getting molecule...\n"
    molecule = Molecule()
    @output "Computing AO Integrals...\n"
    aoint = ConventionalAOIntegrals()

    @output "Calling RHF module...\n"
    refwfn = Fermi.HartreeFock.RHF(molecule, aoint)
    CASCI{T}(refwfn, Alg)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, Alg::ACI; ci = nothing) where T <: AbstractFloat
    @output "Generating Integrals for CAS computation...\n"
    aoint = ConventionalAOIntegrals()
    @output "Transforming Integrals for CAS computation...\n"
    # Read options
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

    s = 1:(frozen+active)

    h = T.(Fermi.Integrals.transform_fock(aoint.T+aoint.V, refwfn.C[:,s], refwfn.C[:,s]))
    V = T.(Fermi.Integrals.transform_eri(aoint.ERI, refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s]))

    aoint = nothing

    #if Fermi.CurrentOptions["cas_use_no"]
    #    eps = refwfn.eps
    #    ndocc = refwfn.ndocc
    #    nvir = refwfn.nvir
    #    den = [i + j - a - b for i in eps[1:ndocc], j in eps[1:ndocc], a in eps[ndocc+1:end], b in eps[ndocc+1:end]]
    #    o = 1:ndocc
    #    v = ndocc+1:ndocc+nvir
    #    VV = T.(Fermi.Integrals.transform_eri(aoint.ERI, refwfn.C[:,o], refwfn.C[:,v], refwfn.C[:,o], refwfn.C[:,v]))
    #    VV = permutedims(VV,(1,3,2,4))
    #    Dn = Array{T}(undef,size(nvir,nvir))
    #    Dd = Array{T}(undef,size(nvir,nvir))
    #    println(size(VV))
    #    @tensor begin 
    #        Dn[a,b] := VV[i,i,c,b]*VV[i,i,c,a]
    #        Dd[a,b] := (den[i,i,c,b]*den[i,i,c,a])
    #    end
    #    D = Dn ./ Dd
    #    refwfn.C[:,v] .= refwfn.C[:,v]*D
    #    display(refwfn.C)
    #    h = T.(Fermi.Integrals.transform_fock(aoint.T+aoint.V, refwfn.C[:,s], refwfn.C[:,s]))
    #    V = T.(Fermi.Integrals.transform_eri(aoint.ERI, refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s]))
    #end


    aoint = nothing
    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg, ci=ci)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::ACI; ci=nothing) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    ttotal = @elapsed begin
    @output "\n    • Computing FCI with the ACI algorithm.\n\n"
    act_range = (frozen+1):(active+frozen)
    σ = Fermi.CurrentOptions["σ"]
    γ = Fermi.CurrentOptions["γ"]
    pe = Fermi.CurrentOptions["aci_print_screen"]

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active

    # Start reference space as HF
    zeroth = repeat('1', frozen)*repeat('1', Int(act_elec/2))
    if ci == nothing
        P = [Determinant(zeroth, zeroth)]
        Pcoef = [1.0]
    else
        P = deepcopy(ci.dets)
        Pcoef = deepcopy(ci.coef)
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
    cflag = true
    while true
        if ite > 20
            break
        end
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
        t = @elapsed Fe = ϵI(F, P, Pcoef, E, h, V)
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
        @output "Updating model space ..." length(M)
        t = @elapsed E, Pcoef, P = update_model_space(M, h, V)
        @output " Model space updated in {:5.5f}\n" t
        ΔE += E
        @output "Model Space Energy           {:15.10f}\n" E + refwfn.molecule.Vnuc
        @output "Energy Change                {:15.10f}\n" ΔE

        if oldP == Set(P) 
            break
        end
        oldP = Set(deepcopy(P))
        Lenny = length(P)
        @output "Coarse graining model space for next iteration\n"
        # Coarse grain
        Cperm = zeros(Int, length(P))
        sortperm!(Cperm, Pcoef, by=i->i^2)
        reverse!(Cperm)

        Pcoef = Pcoef[Cperm]
        P = P[Cperm]

        while true
            if sum(Pcoef[1:end-1].^2) >= 1-γ*σ
                pop!(Pcoef)
                pop!(P)
            else
                break
            end
        end
        @output "Final coarse grained model space size is {}\n" length(P)
        @output repeat("=",50)*"\n"
        ite += 1
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
    E = (E+refwfn.molecule.Vnuc)

    @output "\n • Most important determinants:\n\n"

    for i in 1:(min(20,length(P)))
        @output "{:15.5f}      {}\n" Pcoef[i]  detstring(P[i], frozen+active)
    end
    
    CASCI{T}(refwfn, E, P, Pcoef)
end

@fastmath @inbounds function get_fois(dets::Array{Determinant,1}, Nα::Int, Nβ::Int, act_range::UnitRange{Int64})::Array{Determinant,1}

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
                    newα = (d.α ⊻ (1<<(i-1))) | (1<<(a-1)) 
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
                            newestα = (newα ⊻ (1<<(j-1))) | (1<<(b-1)) 
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
                    newβ = (d.β ⊻ (1<<(i-1))) | (1<<(a-1)) 
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
                            newestβ = (newβ ⊻ (1<<(j-1))) | (1<<(b-1)) 
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
                            newα = (d.α ⊻ (1<<(i-1))) | (1<<(a-1)) 
                            newβ = (d.β ⊻ (1<<(j-1))) | (1<<(b-1)) 
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
    return fois
end

@fastmath @inbounds function ϵI(Fdets, P::Array{Determinant,1}, Pcoef::Array{Float64,1}, Ep::T, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat
    Fe = zeros(length(Fdets))
    N = sum(αlist(P[1]))
    αinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    βinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    @sync for i in eachindex(Fdets)
        Threads.@spawn begin
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

function update_model_space(M::Array{Determinant,1}, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat

    M = complete_set(M)
    H = get_sparse_hamiltonian_matrix(M, h, V, Fermi.CurrentOptions["cas_cutoff"])

    @output "Diagonalizing Hamiltonian...\n"
    decomp, history = partialschur(H, nev=1, tol=10^-12, which=LM())
    λ, ϕ = partialeigen(decomp)
    #λ,ϕ = eigen(Array(H))

    return λ[1], ϕ[:,1], deepcopy(M)
end


function incomp(dets::Array{Determinant,1})

    ming = 0
    for d in dets

        _det = Determinant(d.β, d.α)
        if !(_det in dets)
            ming += 1
        end
    end

    return ming
end

function complete_set(dets::Array{Determinant,1})

    newdets = Determinant[]
    for d in dets
        
        asym = d.α ⊻ d.β
        if asym == 0
            continue
        end

        sym = d.α & d.β

        n = count_ones(asym)
        e = Int(n/2)
        idx = Int[]

        str = repeat("1",e)*repeat("0",e)
        perms = multiset_permutations(str, n)
        
        i = 1
        while i ≤ asym
            if 1<<(i-1) & asym ≠ 0
                push!(idx, i) 
            end
            i += 1
        end

        foo = x->reverse(bitstring(x))[1:7]
        for p in perms
            newα = sym
            newβ = sym
            for (x,i) in zip(p,idx)
                if x == '1'
                    newα = newα | (1<<(i-1))
                elseif x == '0'
                    newβ = newβ | (1<<(i-1))
                end
            end
            push!(newdets, Determinant(newα, newβ))
        end
    end

    return unique(vcat(dets,newdets))
end
