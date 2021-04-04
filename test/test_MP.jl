@reset
@set printstyle none

Econv = [
  -76.330243064527579
  -56.407496860743720
  -231.145778408807132
  -154.579751163193521
  -114.174489244361155
  -279.363387904070862
  -40.427965089111645
  -378.507540337610067
  -228.295786774063146
  -109.338053128796020
]

Edf = [
  -76.33022227893258
  -56.40749546292524
  -231.14577305753187
  -154.57973239300361
  -114.17447470611683
  -279.36344752104367
  -40.42793797440559
  -378.50750865474572
  -228.29574274137781
  -109.33776675040049
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

