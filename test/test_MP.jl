@reset
@set printstyle none

Econv = [
-76.3302430645
-56.4074968607  
-231.1457784088  
-154.5797511632  
-114.1744892444  
-283.2299652413 
-40.4279650891   
-378.5075403376  
-228.295786774110  
-109.3380531288
]

Edf = [
 -76.3302222789  
 -56.4074954629  
 -231.1457730572     
 -154.5797323931
 -114.1744747061
 -283.2299482914
 -40.4279379744  
 -378.5075086542
 -228.2957427415
 -109.3377667504
]

@testset "RMP2" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")
	    
            if molecules[i] == "glycine"
	        continue
            end

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rmp2
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Density Fitted" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvqz-jkfit")
        Fermi.Options.set("rifit", "cc-pvqz-rifit")

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")

            # Skipping these cause there is some problem with Cart x Spherical
            if molecules[i] in ["benzene", "phosphaethene", "glycine"]
                continue
            end

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rmp2
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end

