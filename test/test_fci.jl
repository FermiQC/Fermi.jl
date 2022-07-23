@reset
@set printstyle none

psi4vals = [
    -75.018843443916538,
    -154.015866713653679
]


@testset "FCI" begin

    # WATER STO-3G
    # Read molecule
    path = joinpath(@__DIR__, "xyz/water.xyz")
    mol = open(f->read(f,String), path)
	   
    # Define options
    Fermi.Options.set("molstring", mol)
    Fermi.Options.set("basis", "sto-3g")

    wf = @energy fci
    @test isapprox(wf.energy, psi4vals[1], rtol = 5e-7)

    # ETHANOL FROZEN VIRTUALS
    path = joinpath(@__DIR__, "xyz/ethanol.xyz")
    mol = open(f->read(f,String), path)
	   
    # Define options
    Fermi.Options.set("molstring", mol)
    Fermi.Options.set("basis", "6-31g")
    Fermi.Options.set("drop_vir", 24)

    wf = @energy fci
    @test isapprox(wf.energy, psi4vals[2], rtol = 5e-7)
end