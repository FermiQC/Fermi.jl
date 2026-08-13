@reset
@set printstyle none

# Reference dipole magnitudes (a.u.) from Psi4 1.10, same geometries/bases already
# validated against Fermi's own energies in test_RHF.jl/test_UHF.jl. Psi4 reorients
# molecules by default, so components aren't directly comparable -- the magnitude
# is orientation-independent and is what's checked here.
dipole_tol = 1E-5

@testset "Dipole Moment" begin
    @testset "RHF" begin
        # (molecule index into `molecules`/`basis`, Psi4 reference |μ| in a.u.)
        cases = [(1, 0.8069241612468522),  # water / cc-pvtz
                 (5, 1.084442942581049),   # formaldehyde / 6-31g*
                 (10, 0.0)]                # nitrogen / cc-pcvdz -- zero by symmetry

        for (i, ref) in cases
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rhf
            μ = Fermi.HartreeFock.dipole(wf)
            @test isapprox(√sum(abs2, μ), ref, atol=dipole_tol)
        end
    end

    @testset "UHF" begin
        # OH radical / cc-pvtz
        path = joinpath(@__DIR__, "xyz/oh.xyz")
        mol = open(f->read(f,String), path)

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "cc-pvtz")
        Fermi.Options.set("reference", "uhf")
        Fermi.Options.set("charge", 0)
        Fermi.Options.set("multiplicity", 2)

        wf = @energy uhf
        μ = Fermi.HartreeFock.dipole(wf)
        @test isapprox(√sum(abs2, μ), 0.6973488937715733, atol=dipole_tol)
    end
end

@reset
