@reset
@set printstyle none

Econv = [
-76.330243064527608
-56.407496860743684
-231.145778408807047
-154.579751163193492
-114.167209026284311
-279.363387904071317
-40.427965089111503
-378.507540337610976
-228.295786774063089
-109.338053128796020
]

Edf = [
-76.33022227893399
-56.40749546292774
-231.14576763925319
-154.57973239308316
-114.16719481064435
-279.36344752110256
-40.42793797440178
-378.50753650237135
-228.29574274154695
-109.33776675040524
]

@testset "RMP2" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)
	    
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
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rmp2
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end

