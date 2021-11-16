using GaussianBasis
using Molecules
using TensorOperations

function gradient(wfn::RHF)

    # Following eq. on C.3. Szabo & Ostlund
    atoms = wfn.molecule.atoms
    Natoms = length(atoms)
    bset = BasisSet(wfn.orbitals.basis, atoms)
    nbas = bset.nbas

    ∂E = zeros(Natoms, 3)

    @views Co = wfn.orbitals.C[:,1:wfn.ndocc]

    P = 2.0 * Co * Co'
    Q = 2.0 * Co * diagm(wfn.orbitals.eps[1:wfn.ndocc]) * Co'

    # Preallocate arrays
    ∂H = zeros(nbas, nbas, 3)
    ∂S = zeros(nbas, nbas, 3)
    ∂ERI = zeros(3, nbas, nbas, nbas, nbas)

    for a in eachindex(atoms)

        ∂H .= 0
        ∂S .= 0
        ∂ERI .= 0

        # Nuclear repulsion
        ∂E[a, :] .= Molecules.∇nuclear_repulsion(atoms, a)

        # Store kinetic into H and nuclear into S
        GaussianBasis.∇kinetic!(∂H, bset, a)
        GaussianBasis.∇nuclear!(∂S, bset, a)

        ∂H += ∂S
        ∂S .= 0

        # Now use S for overlap
        GaussianBasis.∇overlap!(∂S, bset, a)

        GaussianBasis.∇ERI_2e4c!(∂ERI, bset, a)

        for q in 1:3
            ∂E[a, q]  = sum(P .* ∂H[:,:,q])
            ∂E[a, q] -= sum(Q .* ∂S[:,:,q])

            ASERI = ∂ERI[q,:,:,:,:] - 0.5*permutedims(∂ERI[q,:,:,:,:], (1,4,3,2))
            @tensoropt X = 0.5*P[μ,ν]*P[λ,σ]*ASERI[μ,ν,σ,λ]
            ∂E[a,q] += X
        end
    end

    return ∂E

end