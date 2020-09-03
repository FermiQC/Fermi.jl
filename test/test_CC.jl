@reset
@testset "Coupled Cluster" begin
    wf = @energy ccsd;
    @test isapprox(wf.CorrelationEnergy, -75.0187095714, rtol=tol)
    wf = @energy ccsd(t);
    @test isapprox(wf.correction, -7.38086348153572e-5, rtol=tol)
    wf = @energy ecCCSD;
    @test isapprox(wf.CorrelationEnergy, -75.0187251416, rtol=tol)
end
