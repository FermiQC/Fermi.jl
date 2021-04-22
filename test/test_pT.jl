@reset
@set printstyle none

Econv = [
    -76.343819598166874,
    -56.427768639264904,
   -231.209805921161546,
   -154.629287842261533,
   -114.197092510628892,
   -279.422940929334800,
    -40.455101412356385,
   -378.537627014158033,
   -228.331709988470351,
   -109.355575399443808
]

Edf = [
]

@testset "CCSD(T)" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")
	    
            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy ccsd
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    #@testset "Density Fitted" begin
    #    Fermi.Options.set("df", true)
    #    Fermi.Options.set("jkfit", "cc-pvqz-jkfit")
    #    Fermi.Options.set("rifit", "cc-pvqz-rifit")

    #    for i = eachindex(molecules)
    #        # Read molecule
    #        mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")

    #        # Skipping these cause there is some problem with Cart x Spherical
    #        if molecules[i] in ["benzene", "phosphaethene"]
    #            continue
    #        end

    #        # Define options
    #        Fermi.Options.set("molstring", mol)
    #        Fermi.Options.set("basis", basis[i])

    #        wf = @energy ccsd
    #        @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
    #    end
    #end
end

