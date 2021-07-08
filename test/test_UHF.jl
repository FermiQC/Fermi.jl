@reset
@set printstyle none

Econv = [-75.4196061687054, -39.5638172309568, -38.9378598408943]

@testset "UHF" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)
        for i = eachindex(uhf_boys)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("scf_max_iter", 500)
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            #Iu = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())

            wf = @energy uhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
end