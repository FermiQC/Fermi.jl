@reset
@set printstyle none
@set lints false

Econv = [
-76.343819598166903
-56.427768639264869
-231.209805921161490
-154.629287842261505
-114.189827180824139
-279.422940929335255
-40.455101412356250
-378.537627014158943
-228.331709988470323
-109.355575399443808
]

Edf = [
-76.34382305394102
-56.42777397519046
-231.20979718974777
-154.62927721611624
-114.18981724346311
-279.42300150553956
-40.45509944258317
-378.53762406241179
-228.33167903449799
-109.35531038066669
]

@testset "CCSD(T)" begin
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

