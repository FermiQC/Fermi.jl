@reset
@set printstyle none

# Psi4, water.xyz, sto-3g, h=0.005
ref_grad = [ 8.3320809487e-03  8.5752642382e-03  8.0196787371e-03
             4.2988605000e-03 -8.6636480397e-03 -8.1043008494e-03
            -1.2629991182e-02  8.7293599531e-05  8.3400498071e-05]

@testset "Tools" begin
    # Finite differences gradients
    @testset "Findif" begin
        Fermi.Options.set("deriv_type", "findif")
        g = @gradient rhf
        Δg = g - ref_grad
        @test √(sum(Δg.^2)/length(g)) < tol
    end
end