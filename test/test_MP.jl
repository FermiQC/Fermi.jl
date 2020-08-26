@reset
@testset "Moller-Plesset" begin
    @set {
          basis cc-pvdz
          scf_alg df
          jkfit cc-pvdz-jkfit
          rifit cc-pvdz-rifit
         }
    wf = @energy rmp2;
    @test isapprox(wf.CorrelationEnergy, -0.2062284914, rtol=tol)
end
