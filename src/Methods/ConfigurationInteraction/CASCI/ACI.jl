using Combinatorics
using SparseArrays
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
    act_range = (frozen+1):active
    σ = Fermi.CurrentOptions["σ"]/1000
    γ = Fermi.CurrentOptions["γ"]/1000

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active

    # Start reference space as HF
    zeroth = repeat('1', frozen)*repeat('1', Int(act_elec/2))
    Pdets = [Determinant(zeroth, zeroth)]
    Pcoef = [1.0]
    Penergy = refwfn.energy - refwfn.molecule.Vnuc
    ΔPenergy = 1.0
    ite = 1

    @output repeat("=",50)*"\n"
    while abs(ΔPenergy) > 10^-7
        if ite > 20
            break
        end
        @output "   Iteration {}\n" ite
        @output "Initial model space size: {}\n" length(Pdets)
        @output "Generating First Order Interacting Space...\n"
        F = get_fois(Pdets, Int(act_elec/2), Int(act_elec/2), act_range)
        @output "FOIS size: {}\n" length(F)
        Fe = ϵI(F, Pdets, Pcoef, Penergy, h, V)
        @output "Screening FOIS for σ = {}\n" σ

        # Sort F by energy contribution
        Fperm = zeros(Int, length(Fe))
        sortperm!(Fperm, Fe, by=abs)
        Fperm = reverse(Fperm)

        Fe = Fe[Fperm]
        F = F[Fperm]

        ϵsum = 0.0
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
        @output "New model space size {}\n" length(Pdets)+length(F)

        # Now we build the total space
        ΔPenergy = -Penergy
        Penergy, Pcoef, Pdets = update_model_space(Pdets, F, h, V)
        ΔPenergy += Penergy
        @output "Model Space Energy           {:15.10f}\n" Penergy + refwfn.molecule.Vnuc
        @output "Energy Change                {:15.10f}\n" ΔPenergy

        @output "Coarse graining model space for next iteration\n"
        # Coarse grain
        Cperm = zeros(Int, length(Pdets))
        sortperm!(Cperm, Pcoef, by=i->i^2)
        Cperm = reverse(Cperm)

        Pcoef = Pcoef[Cperm]
        Pdets = Pdets[Cperm]

        Csum = 0.0
        while true
            Csum += Pcoef[end]^2
            if Csum ≤ 0.001-γ*σ
                pop!(Pcoef)
                pop!(Pdets)
            else
                Csum -= Pcoef[end]^2
                break
            end
        end
        @output repeat("=",50)*"\n"
        ite += 1
    end
end

function get_fois(dets::Array{Determinant,1}, Nα::Int, Nβ::Int, act_range::UnitRange{Int64})

    # Ns must be > 1
    fois = Determinant[]

    αocc = zeros(Int,Nα)
    βocc = zeros(Int,Nβ)
    αuno = zeros(Int,length(act_range)-Nα)
    βuno = zeros(Int,length(act_range)-Nβ)
    for d in dets
        αindex!(d,αocc)
        βindex!(d,βocc)
        αvir!(d, act_range, αuno)
        βvir!(d, act_range, βuno)
        # Get αα -> αα excitations
        for i in αocc
            for a in αuno
                _, _det = annihilate(d,i,'α')
                _, _det = create(_det, a,'α')
                _det in fois || _det in dets ? nothing : push!(fois, _det)
                for j in αocc
                    if j ≥ i
                        break
                    end
                    for b in αuno
                        if b ≥ a
                            break
                        end
                        _, _det = annihilate(d,i,'α')
                        _, _det = annihilate(_det,j,'α')
                        _, _det = create(_det,a,'α')
                        _, _det = create(_det,b,'α')
                        _det in fois || _det in dets ? nothing : push!(fois, _det)
                    end
                end
            end
        end
        # Get ββ -> ββ excitations
        for i in βocc
            for a in βuno
                _, _det = annihilate(d,i,'β')
                _, _det = create(_det, a,'β')
                _det in fois || _det in dets ? nothing : push!(fois, _det)
                for j in βocc
                    if j ≥ i
                        break
                    end
                    for b in βuno
                        if b ≥ a
                            break
                        end
                        _, _det = annihilate(d,i,'β')
                        _, _det = annihilate(_det,j,'β')
                        _, _det = create(_det,a,'β')
                        _, _det = create(_det,b,'β')
                        _det in fois || _det in dets ? nothing : push!(fois, _det)
                    end
                end
            end
        end
        # Get αβ -> αβ excitations
        for i in αocc
            for a in αuno
                for j in βocc
                    for b in βuno
                        _, _det = annihilate(d,i,'α')
                        _, _det = annihilate(_det,j,'β')
                        _, _det = create(_det,a,'α')
                        _, _det = create(_det,b,'β')
                        _det in fois || _det in dets ? nothing : push!(fois, _det)
                    end
                end
            end
        end
    end
    return fois
end

function ϵI(Fdets::Array{Determinant,1}, Pdets::Array{Determinant,1}, Pcoef::Array{Float64,1}, Ep::T, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat
    # Compute energy estimates
    Fe = zeros(length(Fdets))
    N = sum(αlist(Pdets[1]))
    αind = Array{Int64,1}(undef,N)
    βind = Array{Int64,1}(undef,N)
    for i in eachindex(Fdets)
        D1 = Fdets[i]
        αindex!(D1, αind)
        βindex!(D1, βind)
        Ei = Hd0(αind, βind, h, V)
        Δ = Ei - Ep
        Vint = 0.0
        for j in eachindex(Pdets)
            D2 = Pdets[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                continue 
            elseif el == 2
                Vint += Pcoef[j]*Hd2(D1, D2, V, αexc)
            elseif el == 1
                Vint += Pcoef[j]*Hd1(αind, βind, D1, D2, h, V, αexc)
            else
                nothing
            end
        end
        ϵ = Δ/2 - √((Δ^2)/4 + Vint^2)
        Fe[i] = ϵ
    end
    return Fe
end

function update_model_space(P::Array{Determinant,1}, Q::Array{Determinant,1}, h::Array{T,2}, V::Array{T,4}) where T <: AbstractFloat
    M = vcat(P, Q)

    H = get_sparse_hamiltonian_matrix(M, h, V, Fermi.CurrentOptions["cas_cutoff"])

    @output "Diagonalizing Hamiltonian...\n"
    decomp, history = partialschur(H, nev=1, tol=10^-12, which=LM())
    λ, ϕ = partialeigen(decomp)

    return λ[1], ϕ[:,1], M
end
