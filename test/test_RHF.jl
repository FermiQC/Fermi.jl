@reset
@set printstyle none
@testset "RHF" begin
    @testset "Conventional" begin
        wf1 = @energy rhf;
        @test isapprox(wf1.energy, -74.9650028737304410, rtol=tol)
    end
    
    @testset "DF" begin
        @set {
              df true
              basis cc-pvdz
              jkfit cc-pvdz-jkfit
             }
        wf2 = @energy rhf;
        @test isapprox(wf2.energy, -76.02317515357174, rtol=tol)
    end
end
