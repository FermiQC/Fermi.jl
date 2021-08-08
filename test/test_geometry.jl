@testset "Geometry" begin
   
    # Test display
    A1 = Fermi.Geometry.Atom("H", 1, (0.0, 0.0, 0.0))
    x = Fermi.Geometry.string_repr(A1)
    rex = r"Symbol:\s+?H\nAtomic Number:\s+?1\nPosition:\s+?0.0+?\s+?0.0+?\s+?0.0+?"
    @test occursin(rex, x)

    # Test nuclear repulsion
    A2 = Fermi.Geometry.Atom("H", 1, (0.0, 0.0, Fermi.PhysicalConstants.bohr_to_angstrom))
    @test Fermi.Geometry.nuclear_repulsion(A1, A2) == 1.00

    molstringA = """
    C       -2.1315511243      2.2861688237      0.0000000000                 
    H       -1.0615511243      2.2861688237      0.0000000000                 
    H       -2.4882139062      1.4081046164      0.4966839113                 
    H       -2.4882187621      2.2950594327     -1.0087661539                 
    H       -2.4882200570      3.1553408443      0.5120813130"""

    mol = Fermi.Geometry.Molecule(molstring=molstringA, unit="angstrom", charge=0, multiplicity=1)
    @test begin
        length(mol.atoms) == 5 &&
        mol.charge == 0 &&
        mol.multiplicity == 1 &&
        mol.Nα == mol.Nβ == 5 &&
        mol.Vnuc ≈ 13.6865180253
    end

    molstringB = """
    C       -6.7482761629      2.9988214333      0.0000000000                 
    H       -5.0743010691      4.2609493569      0.0000000000                 
    H       -8.4961890559      4.1563877616      0.0000000000 
    """
    mol = Fermi.Geometry.Molecule(molstring=molstringB, unit="bohr", charge=0, multiplicity=3)
    @test begin
        length(mol.atoms) == 3 &&
        mol.charge == 0 &&
        mol.multiplicity == 3 &&
        mol.Nα == 5 &&
        mol.Nβ == 3 &&
        mol.Vnuc ≈  6.0160249492
    end

    # Test wrong charge, mult and unit
    @test_throws Fermi.Options.FermiException Fermi.Geometry.Molecule(molstring=molstringA, charge=10)
    @test_throws Fermi.Options.FermiException Fermi.Geometry.Molecule(molstring=molstringA, charge=0, multiplicity=2)
    @test_throws Fermi.Options.FermiException Fermi.Geometry.Molecule(molstring=molstringA, unit="invalid")

    # Invalid XYZ
    molstringC = """
    H       0.0000000000                 
    H       0.0000000000 
    """
    @test_throws Fermi.Options.FermiException Fermi.Geometry.Molecule(molstring=molstringC)

    molstringD = """
    H       -5.0743010691      4.2609493569      0.0000000000                 
    H       -8.4961890559      HelloThere        0.0000000000
    """
    @test_throws Fermi.Options.FermiException Fermi.Geometry.Molecule(molstring=molstringD)

    # Test print molecule
    mol = Fermi.Geometry.Molecule(molstring=molstringB, unit="bohr", charge=0, multiplicity=3)

    x = Fermi.Geometry.string_repr(mol)

    @test begin
        occursin("Molecule:", x)          &&
        occursin(r"C\s+?[+-]{0,1}\d+?\.\d+?\s+?[+-]{0,1}\d+?\.\d+?\s+[+-]{0,1}\d+?\.\d+?",x) &&
        occursin(r"H\s+?[+-]{0,1}\d+?\.\d+?\s+?[+-]{0,1}\d+?\.\d+?\s+[+-]{0,1}\d+?\.\d+?",x) &&
        occursin(r"Charge:\s+?0", x)      &&
        occursin(r"Multiplicity:\s+?", x) &&
        occursin(r"Nuclear\s+?repulsion:\s+?\d+?\.\d+?", x)
    end
end