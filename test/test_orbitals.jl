@testset "Orbitals" begin
    @set basis cc-pvdz
    orbs = Fermi.Orbitals.GeneralRestrictedOrbitals(rand(5,5), sd_energy=-50.0)
    @test orbs.basis == "cc-pvdz" && orbs.sd_energy == -50.0 && size(orbs.C) == (5,5)
end