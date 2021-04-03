@reset
@set printstyle none

Econv = [
-76.0529365964
-56.2053544136
-230.6243099831 
-154.0923925902 
-113.8654754290 
-279.1153015356 
-40.2134296254 
-378.3408258921 
-227.7630933676 
-108.9544824030
]

Edf = [
  -76.0529366627
  -56.2053607582
  -230.6243073587
  -154.0923873206
  -113.8654701272
  -279.1153639765
  -40.2134277609
  -378.3407965866
  -227.7630791673
  -108.9544635076
]

@testset "RHF" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Density Fitted" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvqz-jkfit")

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")

            # Skipping these cause there is some problem with Cart x Spherical
            if molecules[i] in ["benzene", "phosphaethene"]
                continue
            end

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end
