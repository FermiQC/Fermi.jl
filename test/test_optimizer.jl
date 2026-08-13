@reset
@set printstyle none

using LinearAlgebra

# Reference values from Psi4 1.10 (same geometries/bases used elsewhere in the
# test suite). Tolerances are looser than the tight energy tests elsewhere
# since cross-code optimizer convergence-threshold differences are normal --
# both cases here agreed with Psi4 to ~1e-5/1e-8 level when checked directly.
opt_tol = 1e-4

function distort(mol::Fermi.Molecule, iA::Int, disp)
    atoms = deepcopy(mol.atoms)
    a = atoms[iA]
    atoms[iA] = Fermi.Atom(a.Z, a.mass, collect(a.xyz) .+ disp)
    return Fermi.Molecule(atoms, mol.charge, mol.multiplicity)
end

@testset "Geometry Optimization" begin
    @testset "RHF water (sto-3g)" begin
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)
        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")

        M = Fermi.Molecule()
        dir = normalize(collect(M.atoms[2].xyz - M.atoms[1].xyz))
        distorted = distort(M, 2, 0.15 .* dir)
        distorted = distort(distorted, 3, [0.1, -0.05, 0.05])

        wfn = Fermi.optimize_geometry("rhf", distorted)

        # Psi4 reference: E = -74.96590118039808, r(OH) = 0.9893373391 Ang, angle = 100.0253935814 deg
        @test isapprox(wfn.energy, -74.96590118039808, atol=opt_tol)

        r1 = norm(collect(wfn.molecule.atoms[2].xyz - wfn.molecule.atoms[1].xyz))
        r2 = norm(collect(wfn.molecule.atoms[3].xyz - wfn.molecule.atoms[1].xyz))
        v1 = collect(wfn.molecule.atoms[2].xyz - wfn.molecule.atoms[1].xyz)
        v2 = collect(wfn.molecule.atoms[3].xyz - wfn.molecule.atoms[1].xyz)
        angle = acosd(dot(v1,v2)/(norm(v1)*norm(v2)))

        @test isapprox(r1, 0.9893373391, atol=opt_tol)
        @test isapprox(r2, 0.9893373391, atol=opt_tol)
        @test isapprox(angle, 100.0253935814, atol=5e-3) # angle in degrees, looser bar
    end

    @testset "RHF N2 diatomic (sto-3g)" begin
        path = joinpath(@__DIR__, "xyz/nitrogen.xyz")
        mol = open(f->read(f,String), path)
        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")

        wfn = Fermi.optimize_geometry("rhf")

        # Psi4 reference: E = -107.500654..., r(NN) = 1.133914497154 Ang
        r = norm(collect(wfn.molecule.atoms[2].xyz - wfn.molecule.atoms[1].xyz))
        @test isapprox(r, 1.133914497154, atol=opt_tol)
    end

    @testset "UHF OH radical (sto-3g)" begin
        path = joinpath(@__DIR__, "xyz/oh.xyz")
        mol = open(f->read(f,String), path)
        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")
        Fermi.Options.set("reference", "uhf")
        Fermi.Options.set("charge", 0)
        Fermi.Options.set("multiplicity", 2)

        M = Fermi.Molecule()
        dir = normalize(collect(M.atoms[2].xyz - M.atoms[1].xyz))
        distorted = distort(M, 2, 0.2 .* dir)

        wfn = Fermi.optimize_geometry("uhf", distorted)

        # Psi4 reference: E = -74.3648856845625232, r(OH) = 1.013964623762 Ang
        @test isapprox(wfn.energy, -74.3648856845625232, atol=opt_tol)
        r = norm(collect(wfn.molecule.atoms[2].xyz - wfn.molecule.atoms[1].xyz))
        @test isapprox(r, 1.013964623762, atol=opt_tol)
    end

    @testset "Linear molecule rejection" begin
        # Synthetic 3-atom collinear geometry (like CO2) -- Cartesian gauge-fixing
        # with 3 anchor atoms is not valid for linear molecules (3N-5 DOF, not
        # 3N-6); this must fail loudly rather than silently optimize wrong.
        atoms = [
            Fermi.Atom(6, 12.0, [0.0, 0.0, 0.0]),
            Fermi.Atom(8, 16.0, [1.16, 0.0, 0.0]),
            Fermi.Atom(8, 16.0, [-1.16, 0.0, 0.0]),
        ]
        linear_mol = Fermi.Molecule(atoms, 0, 1)
        @test_throws Fermi.Options.FermiException Fermi.select_anchors(linear_mol.atoms)
    end

    @testset "@optimize macro" begin
        @reset
        @set printstyle none
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)
        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")

        M = Fermi.Molecule()
        distorted = distort(M, 2, 0.1 .* normalize(collect(M.atoms[2].xyz - M.atoms[1].xyz)))

        wfn = @optimize distorted => rhf
        @test isapprox(wfn.energy, -74.96590118039808, atol=opt_tol)
    end
end

@reset
