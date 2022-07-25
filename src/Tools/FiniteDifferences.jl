function create_displacement(mol::Molecule, A::Int, i::Int, h::Real)

    disp = zeros(3)
    disp[i] += h
    new_mol = deepcopy(mol)
    a = mol.atoms[A]
    new_mol.atoms[A] = Atom(a.Z, a.mass, a.xyz + disp)

    return new_mol
end

function apply_gradient(mol, g, d = 0.001)
    Zvals = [A.Z for A = mol.atoms]
    Svals = [A.AtomicSymbol for A = mol.atoms]

    new_atoms = Fermi.Atom[]

    for i = eachindex(Zvals)
        x = mol.atoms[i].xyz[1] - d*g[i, 1]
        y = mol.atoms[i].xyz[2] - d*g[i, 2]
        z = mol.atoms[i].xyz[3] - d*g[i, 3]
        push!(new_atoms, Fermi.Atom(Svals[i], Zvals[i], (x,y,z)))
    end

    return Fermi.Molecule(new_atoms, mol.charge, mol.multiplicity)
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

function opt_test(energy_function; h=0.005, d=0.01)
    @set printstyle none

    scf_Etol  = Options.get("scf_e_conv")
    scf_Dtol  = Options.get("scf_max_rms")
    old_mol = Fermi.Molecule()

    # Central
    wfn = eval(Expr(:call, energy_function, old_mol))
    oldE = wfn.energy
    println(format("Initial Energy: {:15.10f}", oldE))

    ite = 1
    dE = 1
    while abs(dE) > 1e-8
        g = gradient_test(old_mol, energy_function, h)
        new_mol = apply_gradient(old_mol, g, d)
        wfn = eval(Expr(:call, energy_function, new_mol))
        if abs(wfn.e_conv) > scf_Etol || abs(wfn.d_conv) > scf_Dtol
            display(new_mol)
            error("SCF not converged")
        end
        dE = oldE - wfn.energy
        RMS = √(sum(g.^2) / length(g))
        GRMS = geom_rms(old_mol, new_mol)
        println(format("Iter {:3}   Energy: {:15.10f}   ΔE: {:15.10f}   RMS: {:15.10f}   GRMS: {:15.10f}", ite, wfn.energy, dE, RMS, GRMS))
        oldE = wfn.energy
        old_mol = new_mol
        ite += 1
    end

    display(old_mol)
end