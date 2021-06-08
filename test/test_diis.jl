@testset "DIIS" begin
    dm = Fermi.DIIS.DIISManager{Float32, Float32}()
    @test length(dm) == 0
    x = rand(10)
    y = rand(10)
    push!(dm, x, y)
    @test length(dm) == 1
    @test_throws BoundsError Fermi.DIIS.extrapolate(dm)
end
