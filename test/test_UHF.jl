@reset
@set printstyle none

Econv = [-75.4196061687054,  -39.5638172309568,  -38.9378598408943]
Edf =   [-75.41960080317435, -39.56380446089254, -38.93785583198393]
@testset "UHF" begin
    @testset "Conventional Chonky" begin
        Fermi.Options.set("df", false)
        for i = eachindex(uhf_boys)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            I = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())

            wf = @energy I => uhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Conventional SparseERI" begin
        Fermi.Options.set("df", false)
        for i = eachindex(uhf_boys)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            I = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.SparseERI())

            wf = @energy I => uhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end

    @testset "Density Fitting" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvtz-jkfit")
        for i in [1,3]
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])

            wf = @energy uhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end