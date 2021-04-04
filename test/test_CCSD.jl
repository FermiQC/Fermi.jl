@reset
@set printstyle none

Econv = [
  -76.335767822597319
  -56.422272522003759
  -231.188695053088651
  -154.616795317070171
  -114.187588904912914
  -279.415437830677320
  -40.448675124014301
  -378.532023978713198
  -228.312964998766120
  -109.342963378094467
]

Edf = [
 -76.33577146094589
 -56.42227793333655
 -231.18869172128097
 -154.61678476264649
 -114.18757895139244
 -279.41549876119188
 -40.44867306183189
 -378.53199337419994
 -228.31293477894931
 -109.34269895623227
]

@testset "CCSD" begin
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

            wf = @energy ccsd
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end

