using Combinatorics
using SparseArrays
using LinearAlgebra
using ArnoldiMethod

function CASCI{T}(Alg::ACI) where T <: AbstractFloat

    @output "Getting molecule...\n"
    molecule = Molecule()
    @output "Computing AO Integrals...\n"
    aoint = ConventionalAOIntegrals()

    @output "Calling RHF module...\n"
    refwfn = Fermi.HartreeFock.RHF(molecule, aoint)

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
    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::ACI) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the ACI algorithm.\n\n"
    act_range = (frozen+1):(active+frozen)
    σ = Fermi.CurrentOptions["σ"]
    γ = Fermi.CurrentOptions["γ"]
    pe = Fermi.CurrentOptions["aci_print_screen"]

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active

    Hd0s = Dict{Determinant,T}()
    # Start reference space as HF
    zeroth = repeat('1', frozen)*repeat('1', Int(act_elec/2))
    P = [Determinant(zeroth, zeroth)]
    Pcoef = [1.0]
    E = refwfn.energy - refwfn.molecule.Vnuc
    ΔE = 1.0
    ite = 1

    @output repeat("=",50)*"\n"
    Nα = Int(act_elec/2)
    Nβ = Int(act_elec/2)
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
        sleep(1)
        t = @elapsed begin
            Fperm = zeros(Int, length(Fe))
            sortperm!(Fperm, Fe, by=abs)
            #reverse!(Fperm)
            Fperm = reverse(Fperm)
            Fe = Fe[Fperm]
            F = F[Fperm]
        end
        println(Fe)
        @output "Sorted in {:5.5f} s.\n" t
        @output "Filtering F...\n"
        t = @elapsed begin
        ϵsum = 0.0
        esum = sum(abs.(Fe))
        while true
            ϵsum += abs(Fe[end])
            if ϵsum ≤ σ
                pop!(Fe)
                pop!(F)
            else
                ϵsum -= abs(Fe[end])
                break
            end
        end
        end
        @output "size of Q {}\n" length(Fe)
        @output "Secondary space (Q) built in {:5.5f}\n\n" t

        @output "   • M = P ∪ Q\n"
        ΔE = -E
        M = vcat(P, F)
        @output "Model space size: {}\n" length(M)
        t = @elapsed E, Pcoef, P = update_model_space(M, h, V)
        @output "Model space updated in {:5.5f}\n" t
        ΔE += E
        @output "Model Space Energy           {:15.10f}\n" E + refwfn.molecule.Vnuc
        @output "Energy Change                {:15.10f}\n" ΔE

        if abs(ΔE) < σ/10
            break
        end
        @output "Coarse graining model space for next iteration\n"
        # Coarse grain
        Cperm = zeros(Int, length(P))
        sortperm!(Cperm, Pcoef, by=i->i^2)
        Cperm = reverse(Cperm)

        Pcoef = Pcoef[Cperm]
        P = P[Cperm]

        while true
            if sum(Pcoef.^2) >= 1-γ*σ
                c = pop!(Pcoef)
                det = pop!(P)
            else
                break
            end
        end
        @output "Final coarse grained model space size is {}\n" length(P)
        @output repeat("=",50)*"\n"
        ite += 1
    end

    E = (E+refwfn.molecule.Vnuc)
    CASCI{T}(refwfn, E, P, Pcoef)
end

function get_fois(dets::Array{Determinant,1}, Nα::Int, Nβ::Int, act_range::UnitRange{Int64})

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

    #fois = vcat(fois...)
    #fois = fois[lst]
    println("Pre size $(length(fois))")
    fois = filter((x)->x != Determinant(0,0), fois)
    println("post1 size $(length(fois))")
    fois = unique(fois)
    println("post2 size $(length(fois))")
    fois = filter((x)->!(x in dets), fois)
    println("post size $(length(fois))")
    return fois
end

function ϵI(Fdets::Array{Determinant,1}, P::Array{Determinant,1}, Pcoef::Array{Float64,1}, Ep::T, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat
    Fe = zeros(length(Fdets))
    N = sum(αlist(P[1]))
    αinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    βinds = [Array{Int64,1}(undef,N) for i=1:Threads.nthreads()]
    tasks = Array{Task}(undef,length(Fe))
    for i in eachindex(Fdets)
        tasks[i] = Threads.@spawn begin
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
                #Vint += Hd2(D1, D2, V, αexc)
            elseif el == 1
                Vint += Pcoef[j]*Hd1(αind, βind, D1, D2, h, V, αexc)
                #Vint += Hd1(αind, βind, D1, D2, h, V, αexc)
            end
        end
        test = zeros(2,2)
        test[1,1] = Ep
        test[2,2] = Ei
        test[1,2] = Vint
        test[2,1] = Vint
        ϵ = Δ/2 - √((Δ^2)/4 + Vint^2)
        val,vec = eigen(Symmetric(test))
        Fe[i] = val[1] - Ep
        end #Threads.@spawn
    end
    wait.(tasks)
    return Fe
end

function update_model_space(M::Array{Determinant,1}, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat

    H = get_sparse_hamiltonian_matrix(M, h, V, Fermi.CurrentOptions["cas_cutoff"])

    @output "Diagonalizing Hamiltonian...\n"
    decomp, history = partialschur(H, nev=1, tol=10^-12, which=LM())
    λ, ϕ = partialeigen(decomp)

    return λ[1], ϕ[:,1], deepcopy(M)
end
