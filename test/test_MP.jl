@reset
@set printstyle none

Econv = [
  -76.3302176134 
  -56.4074869825
  -231.1457760762 
  -154.5797414699 
  -114.1744801956 
  -279.3633858666 
  -40.4279410929 
  -378.5075381588 
  -228.2957609516 
  -109.3377924036 
]

Edf = [
  -76.3302222789
  -56.4074954629
  -231.1457730570
  -154.5797323931
  -114.1744747061
  -279.3634475209
  -40.4279379744
  -378.5075086547
  -228.2957427415
  -109.3377667504
]

@testset "RMP2" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            mol = open(f->read(f,String), "xyz/"*molecules[i]*".xyz")

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
            if molecules[i] in ["benzene", "phosphaethene"]
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

