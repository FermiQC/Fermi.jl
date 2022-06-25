@reset
@set printstyle none

# Psi4, water.xyz, sto-3g, h=0.005
ref_grad = [-0.0340838396 -0.0352244219 -0.0329517863
             0.0487762269  0.0012213568  0.0011425074
            -0.0146915349  0.0340026633  0.0318087599]

@testset "Tools" begin
    # Finite differences gradients
    @testset "Findif" begin
        g = Fermi.gradient_findif(Fermi.HartreeFock.RHF)
        Δg = g - ref_grad
        @test √(sum(Δg.^2)/length(g)) < tol
    end
end