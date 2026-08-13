function create_displacement(mol::Molecule, A::Int, i::Int, h::Real)

    disp = zeros(3)
    disp[i] += h
    new_mol = deepcopy(mol)
    a = mol.atoms[A]
    new_mol.atoms[A] = Atom(a.Z, a.mass, a.xyz + disp)

    return new_mol
end

function geom_rms(mol1, mol2)
    out = 0.0
    N = length(mol1.atoms)
    for a in 1:N
        out += sum((mol1.atoms[a].xyz .- mol2.atoms[a].xyz).^2)
    end
    return √(out / 3*N)
end

function findif_intgrad(X::String, mol, A, i, h=0.005)
    mol_disp = create_displacement(mol, A, i, h)
    I = Fermi.Integrals.IntegralHelper(molecule=mol_disp)
    Xplus = I[X]
    mol_disp = create_displacement(mol, A, i, -h)
    I = Fermi.Integrals.IntegralHelper(molecule=mol_disp)
    Xminus = I[X]
    g = (Xplus - Xminus) ./ (2*h)
    return g * PhysicalConstants.bohr_to_angstrom
end

function gradient_findif(energy_function)
    h = Options.get("findif_disp_size")
    gradient_findif(energy_function, Molecule(), h)
end

function gradient_findif(energy_function, mol::Molecule, h=0.005)
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

    g = (Eplus - Eminus) ./ (2*h)

    return g * PhysicalConstants.bohr_to_angstrom
end