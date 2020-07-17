using Combinatorics
using SparseArrays
using ArnoldiMethod

function CASCI{T}(Alg::SparseHamiltonian) where T <: AbstractFloat

    @output "Getting molecule...\n"
    molecule = Molecule()
    @output "Computing AO Integrals...\n"
    aoint = ConventionalAOIntegrals()

    @output "Calling RHF module..."
    refwfn = Fermi.HartreeFock.RHF(molecule, aoint)

    @output "Transforming Integrals..."
    # Read options
    frozen = Fermi.CurrentOptions["cas_frozen"]

    nmo = refwfn.ndocc + refwfn.nvir

    act_elec = 2*(refwfn.ndocc - frozen)

    if act_elec ≤ 0
        error("Invalid number of frozen orbitals ($frozen) for $(2*refwfn.ndocc) electrons.")
    end

    # Active = -1 means FCI, with frozen
    if Fermi.CurrentOptions["cas_active"] == -1
        active = nmo - frozen
    else
        active = Fermi.CurrentOptions["cas_active"]
    end

    if active ≤ act_elec/2
        error("Number of active orbitals ($active) too small for $(act_elec) active electrons")
    end

    s = 1:(frozen+active)
    h = T.(Fermi.Integrals.transform_fock(aoint.T+aoint.V, refwfn.C[:,s], refwfn.C[:,s]))
    V = T.(Fermi.Integrals.transform_eri(aoint.ERI, refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s], refwfn.C[:,s]))

    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::SparseHamiltonian) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the SparseMatrix algorithm.\n\n"


    nroot = Fermi.CurrentOptions["cas_nroot"]

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active
    
    dets = get_determinants(act_elec, active, frozen)
    Ndets = length(dets)
    @output "Number of Determinants: {:10d}\n" Ndets

    @output "Building Sparse Hamiltonian...\n"

    @time begin
        H = get_sparse_hamiltonian_matrix(dets, h, V, Fermi.CurrentOptions["cas_cutoff"])
    end
    @output "Hamiltonian Matrix size: {:10.6f} Mb\n" Base.summarysize(H)/10^6

    @output "Diagonalizing Hamiltonian for {:3d} eigenvalues...\n" nroot
    @time begin
        decomp, history = partialschur(H, nev=nroot, tol=10^-9, which=LM())
        λ, ϕ = partialeigen(decomp)
    end
    @output "\n Final FCI Energy: {:15.10f}\n" λ[1]+refwfn.molecule.Vnuc

    return CASCI{T}(λ[1]+T(refwfn.molecule.Vnuc), dets, ϕ[:,1])

end

function get_determinants(Ne::Int, No::Int, nfrozen::Int)

    Nae = Int(Ne/2)
    occ_string = repeat('1', nfrozen)

    zeroth = repeat('1', Nae)*repeat('0', No-Nae)

    perms = multiset_permutations(zeroth, length(zeroth))

    dets = Determinant[]
    for αstring in perms
        for βstring in perms
            α = occ_string*join(αstring)
            β = occ_string*join(βstring)
            _det = Determinant(α, β)
            push!(dets, _det)
        end
    end

    # Determinant list is sorted by its excitation level w.r.t the first determinant (normally HF)
    sort!(dets, by=d->excitation_level(dets[1], d))

    return dets
end

function get_sparse_hamiltonian_matrix(dets::Array{Determinant,1}, h::Array{T,2}, V::Array{T,4}, tol::Float64) where T <: AbstractFloat

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))

    αind = Array{Int64,1}(undef,Nα)
    βind = Array{Int64,1}(undef,Nβ)
    vals = T[]
    ivals = Int64[]
    jvals = Int64[]
    elem = 0.0
    tHd0 = 0
    tHd1 = 0
    tHd2 = 0
    t = 0

    for i in 1:Ndets
        D1 = dets[i]
        αind .= αindex(D1, Nα)
        βind .= βindex(D1, Nβ)
        for j in i:Ndets
            D2 = dets[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                nothing
            elseif el == 2
                t = @elapsed elem = Hd2(D1, D2, V, αexc)
                if elem > tol || -elem > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
                tHd2 += t
            elseif el == 1
                t = @elapsed elem = Hd1(αind, βind, D1, D2, h, V, αexc)
                if elem > tol || -elem > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
                tHd1 += t
            else
                t = @elapsed elem = Hd0(αind, βind, h, V)
                if elem > tol || -elem > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
                tHd0 += t
            end
        end
    end

    @output "\n Total Time spent in Hd0: {:10.5f}" tHd0
    @output "\n Total Time spent in Hd1: {:10.5f}" tHd1
    @output "\n Total Time spent in Hd2: {:10.5f}\n" tHd2

    return Symmetric(sparse(ivals, jvals, vals))
end
