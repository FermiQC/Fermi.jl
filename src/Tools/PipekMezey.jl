using LinearAlgebra

function mulliken_charges(C0, ints)


    S = ints["S"]
    q = diag(C0*transpose(C0)*S)

    ## THis might change soon
    basis = ints.orbitals.basisset
    natoms = basis.natoms
    basis_per_atom = zeros(Int, natoms)
    charges = zeros(natoms)
    for a in 1:natoms
        n = 0
        for b in basis[a]
            n += 2*b.l + 1
        end
        basis_per_atom[a] = n
    end

    off = 0
    for a in 1:natoms
        n = basis_per_atom[a]
        charges[a] = sum(q[(off+1):(off+n)])
        off += n
    end

    return charges
end