using GaussianBasis
using Molecules
using TensorOperations

function RHFgrad(x...)
    RHFgrad(Molecule(), x...)
end

function RHFgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        RHFgrad(RHF(mol), x...)
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.HartreeFock.RHF, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

function RHFgrad(wfn::RHF, eri_type::Fermi.Integrals.Chonky)

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
    ∂ERI = zeros(nbas, nbas, nbas, nbas, 3)

    # Intermediate auxiliary arrays
    AUX_ERI = zeros(nbas, nbas, nbas, nbas)
    AUX = zeros(nbas, nbas)

    for a in eachindex(atoms)

        ∂H .= 0
        ∂S .= 0
        ∂ERI .= 0

        # Nuclear repulsion
        ∂E[a, :] .= Molecules.∇nuclear_repulsion(atoms, a)

        # Store kinetic into H and nuclear into S
        GaussianBasis.∇kinetic!(∂H, bset, a)
        GaussianBasis.∇nuclear!(∂S, bset, a)

        ∂H .+= ∂S
        ∂S .= 0

        # Now use S for overlap
        GaussianBasis.∇overlap!(∂S, bset, a)

        GaussianBasis.∇ERI_2e4c!(∂ERI, bset, a)

        for q in 1:3
            @views vH = ∂H[:,:,q]
            AUX .= P .* vH
            ∂E[a, q] += sum(AUX)

            @views vS = ∂S[:,:,q]
            AUX .= Q .* vS
            ∂E[a, q] -= sum(AUX)

            @views ERIq = ∂ERI[:,:,:,:,q]
            permutedims!(AUX_ERI, ERIq, (1,4,3,2))
            AUX_ERI .*= -0.5
            AUX_ERI .+= ERIq

            @tensoropt X = 0.5*P[μ,ν]*P[λ,σ]*AUX_ERI[μ,ν,σ,λ]
            ∂E[a,q] += X
        end
    end

    return ∂E
end

function RHFgrad(wfn::RHF)

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

    # Intermediate auxiliary arrays
    AUX = zeros(nbas, nbas)

    for a in eachindex(atoms)

        ∂H .= 0
        ∂S .= 0

        # Nuclear repulsion
        ∂E[a, :] .= Molecules.∇nuclear_repulsion(atoms, a)

        # Store kinetic into H and nuclear into S
        GaussianBasis.∇kinetic!(∂H, bset, a)
        GaussianBasis.∇nuclear!(∂S, bset, a)

        ∂H .+= ∂S
        ∂S .= 0

        # Now use S for overlap
        GaussianBasis.∇overlap!(∂S, bset, a)

        idx, xyz... = GaussianBasis.∇sparseERI_2e4c(bset, a)

        for q in 1:3
            @views vH = ∂H[:,:,q]
            AUX .= P .* vH
            ∂E[a, q] += sum(AUX)

            @views vS = ∂S[:,:,q]
            AUX .= Q .* vS
            ∂E[a, q] -= sum(AUX)

            ∇q = xyz[q]

            for i in eachindex(idx)
                μ,ν,λ,σ = idx[i] .+ 1
                ∇k = ∇q[i]

                if abs(∇k) < 1e-12
                    continue
                end

                # perms is 1 if all indexes are the same (which must not happen)
                # because it implies that the integral is on only one center e.g. (11|11)
                # for which the derivative must be zero
                # In turn, perms will be 2, 4, or 8 depending on the number of permutations
                # allowed. e.g. (11|23) yields 4: (σ != λ) and μν != λσ
                perms = 1 << (μ !== ν) << (σ !== λ) << (Set((μ,ν)) != Set((σ,λ)))

                if perms === 8
                    ∂E[a,q] += (4.0*P[μ,ν]*P[λ,σ] - P[μ,σ]*P[λ,ν] - P[μ,λ]*P[ν,σ])*∇k
                elseif perms === 4
                    ∂E[a,q] += (2.0*P[μ,ν]*P[λ,σ] - 0.5*(P[μ,σ]*P[λ,ν] + P[μ,λ]*P[ν,σ]))*∇k
                else
                    ∂E[a,q] += (P[μ,ν]*P[λ,σ] - 0.25*(P[μ,σ]*P[λ,ν] + P[μ,λ]*P[ν,σ]))*∇k
                end
            end
        end
    end

    return ∂E
end