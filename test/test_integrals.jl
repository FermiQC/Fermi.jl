@testset "Integrals" begin
    @reset

    @set {
        jkfit auto
        rifit auto
        basis cc-pvtz
    }

    jkfit = Fermi.Integrals.JKFIT()
    rifit = Fermi.Integrals.RIFIT()

    # Test auto-complete
    @test jkfit.basisset.basis_name == "cc-pvtz-jkfit"
    @test rifit.basisset.basis_name == "cc-pvtz-rifit"
end