@reset
@testset "CI" begin
    wf = @energy ci
    @test isapprox(wf.energy,-75.01775104333774,rtol=tol)
    @set ci_alg sparse
    wf = @energy ci
    @test isapprox(wf.energy,-75.0188434735,rtol=tol)
end
