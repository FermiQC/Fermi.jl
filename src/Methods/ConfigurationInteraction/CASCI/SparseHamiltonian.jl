using TensorOperations
using Combinatorics
using SparseArrays
using ArnoldiMethod

function CASCI{T}(Alg::SparseHamiltonian) where T <: AbstractFloat

    @output "Getting molecule...\n"
    molecule = Molecule()
    @output "Computing AO Integrals...\n"
    #aoint = ConventionalAOIntegrals()

    @output "Calling RHF module...\n"
    refwfn = Fermi.HartreeFock.RHF(molecule)
    ints = refwfn.ints

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
    h = T.(Fermi.Integrals.transform_fock(ints["T"] + ints["V"], ints.orbs["FU"][s], ints.orbs["FU"][s]))
    V = T.(Fermi.Integrals.transform_eri(ints["μ"], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s]))

    aoint = nothing
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
    
    @output "\nGenerating determinants ..."
    t = @elapsed begin
        dets = get_determinants(act_elec, active, frozen)
        Ndets = length(dets)
    end
    
    @output " done in {} seconds.\n" t
    @output "\nNumber of Determinants: {:10d}\n" Ndets

    #@output "\nBuilding Sparse Hamiltonian..."

    #t = @elapsed H = get_sparse_hamiltonian_matrix(dets, h, V, Fermi.CurrentOptions["cas_cutoff"])
    #@output " done in {:5.5f} seconds.\n" t

    #@output "Hamiltonian Matrix size: {:10.3f} Mb\n" Base.summarysize(H)/10^6

    #@output "Diagonalizing Hamiltonian for {:3d} eigenvalues..." nroot
    #t = @elapsed begin
    #    decomp, history = partialschur(H, nev=nroot, tol=10^-12, which=LM())
    #    λ, ϕ = partialeigen(decomp)
    #end
    #@output " done in {:5.5f} seconds.\n" t
    #@output "\n Final FCI Energy: {:15.10f}\n" λ[1]+refwfn.molecule.Vnuc


    @output "\nBuilding Strings..."

    t = @elapsed begin
        sts,tree = get_αstrings(act_elec, active)
    end
    @output " done in {:5.5f} seconds.\n" t
    @output "\nBuilding Hamiltonian..."
    t = @elapsed begin
        H = sparse(get_H_fromstrings(sts, tree, h, V, frozen))
        droptol!(H,10^-12) 
    end
    @output " done in {:5.5f} seconds.\n" t

    @output "Hamiltonian Matrix size: {:10.3f} Mb\n" Base.summarysize(H)/10^6

    @output "Diagonalizing Hamiltonian for {:3d} eigenvalues..." nroot
    t = @elapsed begin
        decomp, history = partialschur(H, nev=nroot, tol=10^-12, which=LM())
        λ, ϕ = partialeigen(decomp)
    end
    @output " done in {:5.5f} seconds.\n" t
    @output "\n Final FCI Energy: {:15.10f}\n" λ[1]+refwfn.molecule.Vnuc

    # Sort dets by importance
    C = ϕ[:,1]
    sp = sortperm(abs.(ϕ[:,1]), rev=true)
    C = ϕ[:,1][sp]
    dets = dets[sp]

    @output "\n • Most important determinants:\n\n"

    for i in 1:(min(20,length(C)))
        @output "{:15.5f}      {}\n" C[i]  detstring(dets[i], frozen+active)
    end
    @output "\n"
    return CASCI{T}(refwfn, λ[1]+T(refwfn.molecule.Vnuc), dets, C)
end

function get_determinants(Ne::Int, No::Int, nfrozen::Int)

    Nae = Int(Ne/2)
    occ_string = repeat('1', nfrozen)

    zeroth = repeat('1', Nae)*repeat('0', No-Nae)

    perms = multiset_permutations(zeroth, length(zeroth))

    dets = [Determinant[] for i = 1:Threads.nthreads()]
    @sync for αstring in perms
        Threads.@spawn begin
        for βstring in perms
            α = occ_string*join(αstring)
            β = occ_string*join(βstring)
            _det = Determinant(α, β)
            push!(dets[Threads.threadid()], _det)
        end
        end #Threads.@spawn
    end

    dets = vcat(dets...)
    # Determinant list is sorted by its excitation level w.r.t the first determinant (normally HF)
    #sort!(dets, by=d->excitation_level(dets[1], d))

    return dets
end

function get_sparse_hamiltonian_matrix(dets::Array{Determinant,1}, h::Array{T,2}, V::Array{T,4}, tol::Float64) where T <: AbstractFloat

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))

    αind = [Array{Int64,1}(undef,Nα) for i = 1:Threads.nthreads()]
    βind = [Array{Int64,1}(undef,Nβ) for i = 1:Threads.nthreads()]
    vals = [T[] for i = 1:Threads.nthreads()]
    ivals = [Int64[] for i = 1:Threads.nthreads()]
    jvals = [Int64[] for i = 1:Threads.nthreads()]

    @sync for i in 1:Ndets
        Threads.@spawn begin
        D1 = dets[i]
        αindex!(D1, αind[Threads.threadid()])
        βindex!(D1, βind[Threads.threadid()])
        elem = 0.0
        for j in i:Ndets
            D2 = dets[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                continue 
            elseif el == 2
                elem = Hd2(D1, D2, V, αexc)
            elseif el == 1
                elem = Hd1(αind[Threads.threadid()], βind[Threads.threadid()], D1, D2, h, V, αexc)
            else
                elem = Hd0(αind[Threads.threadid()], βind[Threads.threadid()], h, V)
            end
            if abs(elem) > tol
                push!(vals[Threads.threadid()], elem)
                push!(ivals[Threads.threadid()], i)
                push!(jvals[Threads.threadid()], j)
            end
        end
    end #Threads.@spawn
    end

    ivals = vcat(ivals...)
    jvals = vcat(jvals...)
    vals  = vcat(vals...)
    return Symmetric(sparse(ivals, jvals, vals))
end

function get_1p_coupling_coefficients(dets::Array{Determinant,1}, nmo::Int)

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))
    αind = Array{Int64,1}(undef,Nα)
    βind = Array{Int64,1}(undef,Nβ)    
    γ = zeros(nmo, nmo, Ndets, Ndets)

    for nI in 1:Ndets
        I = dets[nI]
        for nJ in 1:Ndets
            J = dets[nJ]

            if I.α == J.α && I.β == J.β
                αindex!(I,αind) 
                βindex!(I,βind) 

                for i=αind
                    γ[i,i,nI,nI] += 1
                end

                for i=βind
                    γ[i,i,nI,nI] += 1
                end
            end

            αexc = αexcitation_level(I,J)
            if αexc > 1
                continue
            end

            βexc = βexcitation_level(I,J)
            if αexc + βexc > 1
                continue
            end

            if αexc == 1
                i, = αexclusive(I,J)
                j, = αexclusive(J,I)
                p = phase(I,J)
                γ[i,j,nI,nJ] = p

            elseif βexc == 1
                i, = βexclusive(I,J)
                j, = βexclusive(J,I)
                p = phase(I,J)
                γ[i,j,nI,nJ] = p
            end
        end
    end

    return γ
end

function build_H_fullγ(dets, h, V)
    nmo = size(h,1)
    γ = get_1p_coupling_coefficients(dets, nmo)
    δ = [i==j ? 1 : 0 for i = 1:nmo, j = 1:nmo]
    @tensoropt H[I,J] := γ[i,j,I,J]*h[i,j] - 0.5*γ[i,l,I,J]*δ[j,k]*V[i,j,k,l]
    println("One piece")
    display(H)
    @tensoropt H[I,J] += + 0.5*γ[i,j,I,K]*γ[k,l,K,J]*V[i,j,k,l]
    println("Final")
    return H
end
