function create_displacement(mol, A::Int, i::Int, h)

    new_atoms = [mol.atoms...]
    atomA = mol.atoms[A]
    new_xyz = [atomA.xyz...] .+ h
    new_atoms[A] = Fermi.Geometry.Atom(atomA.AtomicSymbol, atomA.Z, (new_xyz[1], new_xyz[2], new_xyz[3]))

    return Fermi.Geometry.Molecule(new_atoms, mol.charge, mol.multiplicity)
end

function gradient_test(energy_function, h=0.005)
    mol = Fermi.Geometry.Molecule()
    N = length(mol.atoms)
    Eplus = zeros(N,3)
    Eminus = zeros(N,3)

    # Plus
    for A = 1:N, i = 1:3
        mol_disp = create_displacement(mol, A, i, h)
        wfn = eval(Expr(:call, energy_function, mol_disp))
        Eplus[A,i] = wfn.energy
    end
    # Minus
    for A = 1:N, i = 1:3
        mol_disp = create_displacement(mol, A, i, -h)
        wfn = eval(Expr(:call, energy_function, mol_disp))
        Eminus[A,i] = wfn.energy
    end

    display(Eplus)
    display(Eminus)
    g = (Eplus - Eminus) ./ (2*h)

    return g
end