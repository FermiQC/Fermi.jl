module FCI
using Combinatorics
using LinearAlgebra
using SparseArrays
using Arpack
using TensorOperations
using Fermi
using Fermi.Output
using Fermi.Wavefunction
using Fermi.ConfigurationInteraction
using Fermi.ConfigurationInteraction.DetOperations
using Fermi.ConfigurationInteraction.MatrixElement

export do_fci

function do_fci(wfn::Wfn; kwargs...)

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the FCI module.\n\n"

    # Process options
    for arg in keys(Fermi.ConfigurationInteraction.defaults)
        if arg in keys(kwargs)
            @eval $arg = $(kwargs[arg])
        else
            @eval $arg = $(Fermi.ConfigurationInteraction.defaults[arg])
        end
    end

    nmo = wfn.nmo

    if wfn.nalpha ≠ wfn.nbeta
        error("Open-shell case not implemented yet.")
    end

    act_elec = wfn.nalpha + wfn.nbeta - 2*frozen

    if act_elec ≤ 0
        error("Invalid number of frozen orbitals ($frozen) for $(nα+nβ) electrons.")
    end

    # Active = -1 means FCI, with frozen
    if active == -1
        @eval active = $nmo - $frozen
    end

    if active ≤ act_elec/2
        error("Number of active orbitals ($active) too small for $(act_elec) activ electrons")
    end

    @output "Transforming integrals...\n"

    # Integral transformation
    C = wfn.Ca
    gao = wfn.ao_eri
    hao = wfn.hao
    @tensoropt begin
        V[i,a,j,b] := C[σ,b]*C[λ,j]*C[ν,a]*C[μ,i]*gao[μ,ν,λ,σ]
        h[p,q] := C[μ,p]*C[ν,q]*hao[μ,ν] 
    end

    @output "Done.\n"

    @output "NMO: {:3d}\n" wfn.nmo

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active
    
    dets = get_determinants(act_elec, active, nmo, frozen)
    Ndets = length(dets)
    @output "Number of Determinants: {:10d}\n" Ndets

    @output "Building Sparse Hamiltonian...\n"

    @time begin
        H = get_sparse_hamiltonian_matrix(dets, h, V, min_matrix_elem)
    end
    @output "Hamiltonian Matrix size: {:10.6f} Mb\n" Base.summarysize(H)/10^6

    @output "Diagonalizing Hamiltonian for {:3d} eigenvalues...\n" nroot
    @time λ, Φ = eigs(H, nev=nroot)
    @output "\n Final FCI Energy: {:15.10f}\n" λ[1]+wfn.vnuc
    #@output "Time: {:10.5f}\n" t

end

function get_determinants(Ne::Int, No::Int, nmo::Int, nfrozen::Int)

    Nae = Int(Ne/2)
    occ_string = repeat('1', nfrozen)
    vir_string = repeat('0', nmo-nfrozen-Nae)

    zeroth = repeat('1', Nae)*repeat('0', No-Nae)

    perms = multiset_permutations(zeroth, length(zeroth))

    dets = Determinant[]
    for αstring in perms
        for βstring in perms
            α = occ_string*join(αstring)*vir_string
            β = occ_string*join(βstring)*vir_string
            _det = Determinant(α, β)
            push!(dets, _det)
        end
    end

    # Determinant list is sorted by its excitation level w.r.t the first determinant (normally HF)
    sort!(dets, by=d->excitation_level(dets[1], d))

    return dets
end


function get_dense_hamiltonian_matrix(dets::Array{Determinant,1}, H::AbstractArray{Float64,2}, h::Array{Float64,2}, V::Array{Float64,4})

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))

    αind = Array{Int64,1}(undef,Nα)
    βind = Array{Int64,1}(undef,Nβ)
    for i in 1:Ndets
        D1 = dets[i]
        αind .= αindex(D1, Nα)
        βind .= βindex(D1, Nβ)
        for j in i:Ndets
            D2 = dets[j]
            el = excitation_level(D1,D2)
            if el > 2
                nothing
            elseif el == 2
                H[i,j] = Hd2(D1, D2, V)
            elseif el == 1
                H[i,j] = Hd1(αind, βind, D1, D2, h, V)
            else
                H[i,j] = Hd0(αind, βind, h, V)
            end
        end
    end

    return Symmetric(H)
end

function get_sparse_hamiltonian_matrix(dets::Array{Determinant,1}, h::Array{Float64,2}, V::Array{Float64,4}, tol::Float64)

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))

    αind = Array{Int64,1}(undef,Nα)
    βind = Array{Int64,1}(undef,Nβ)
    vals = Float64[]
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

function print_matrix(M)

    m,n = size(M)

    for i in 1:m
        println(M[i,:])
    end
end

function sparsity(M::AbstractArray)

    return count(i->abs(i)>1e-10, M)/length(M)

end

end #module
