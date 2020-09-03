@reset
@testset "Coupled Cluster" begin
    wf = @energy CCSD;
    @test isapprox(wf.CorrelationEnergy, -75.0187095714, rtol=tol)
    wf = @energy CCSD(T);
    @test isapprox(wf.correction, -7.38086348153572e-5, rtol=tol)
    wf = @energy ecCCSD;
    @test isapprox(wf.CorrelationEnergy, -75.0187251416, rtol=tol)
    wf = @energy ecCCSD(T);
    @test isapprox(wf.correction, -7.382789470598815e-5, rtol=tol)
    @set ci_alg sparse
    wf = @energy ecCCSD;
    @test isapprox(wf.CorrelationEnergy, -75.0188434735, rtol=tol)
end
