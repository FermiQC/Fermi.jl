@reset
@set printstyle none
@set lints false

Econv = [
-76.05293659641389
-56.20535441552476
-230.62430998333460
-154.09239259145326
-113.86482635166755
-279.11530153592219
-40.21342962572409
-378.34082589236215
-227.76309337209935
-108.95448240313179
]

Edf = [
-76.05293666271527
-56.20536076021406
-230.62430445834565
-154.09238732189544
-113.86482178236324
-279.11536397687240
-40.21342776121311
-378.34082725476162
-227.76307917196510
-108.95446350768205
]

@testset "RHF" begin
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

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
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
        end
    end
end
