@testset "Physical Constants" begin
    x = rand(keys(Fermi.PhysicalConstants.atom_num))
    @test Fermi.PhysicalConstants.atomic_number(x) == Fermi.PhysicalConstants.atom_num[x]
    @test_throws KeyError Fermi.PhysicalConstants.atomic_number("Invalid") 

    @test Fermi.PhysicalConstants.num_core_electrons("Ca") == 18
end
