@reset
@set printstyle none

Econv = [
-76.335767822597347
-56.422272522003723
-231.188695053088594
-154.616795317070142
-114.180708994251702
-279.415437830677774
-40.448675124014166
-378.532023978714108
-228.312964998766091
-109.342963378094467
]

Edf = [
-76.33577146094731
-56.42227793333907
-231.18868648219075
-154.61678476272598
-114.18069925474244
-279.41549876125123
-40.44867306182801
-378.53202120936407
-228.31293477911845
-109.34269895623703
]

@testset "CCSD" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)
	    
            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy ccsd
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Density Fitted" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvqz-jkfit")
        Fermi.Options.set("rifit", "cc-pvqz-rifit")

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)


            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy ccsd
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end

