@reset
@set printstyle none

Econv = [
-76.05293659641389
-56.20535441552476
-230.62430998333460
-154.09239259145326
-113.86482635166755
-279.11530153592219
-40.21342962572409
-378.34082589236215
-227.76309337209935
-108.95448240313179
]

Edf = [
-76.05293666271527
-56.20536076021406
-230.62430445834565
-154.09238732189544
-113.86482178236324
-279.11536397687240
-40.21342776121311
-378.34082725476162
-227.76307917196510
-108.95446350768205
]

ref_grads = Dict{String, Matrix}(
    "water" => [-0.034110486707    -0.035228598159    -0.032954172280
                 0.048800781314     0.001221100923     0.001142264303
                -0.014690294605     0.034007497234     0.031811907977],
    "ammonia" => [0.000000000000     0.000000427189     0.010891836599
                  0.000000000000     0.004413532323    -0.003630495352
                  0.003822520996    -0.002206979756    -0.003630670623
                 -0.003822520996    -0.002206979756    -0.003630670623],
    "benzene" => [-0.000000000000     0.000251826467     0.000000000000
                   0.000217602294     0.000124532780     0.000000000000
                   0.000217602294    -0.000124532780     0.000000000000
                   0.000000000000    -0.000251826467     0.000000000000
                  -0.000217602294    -0.000124532780     0.000000000000
                  -0.000217602294     0.000124532780     0.000000000000
                  -0.000000000000     0.003518154246     0.000000000000
                   0.003047034907     0.001759262144     0.000000000000
                   0.003047034907    -0.001759262144     0.000000000000
                   0.000000000000    -0.003518154246     0.000000000000
                  -0.003047034907    -0.001759262144     0.000000000000
                  -0.003047034907     0.001759262144     0.000000000000],
    "ethanol" => [ 0.000268478447    -0.000722077624     0.000000000000
                   0.004229353933     0.007442087299     0.000000000000
                   0.004195039086    -0.013841867069     0.000000000000
                  -0.007562462384     0.007854197802     0.000000000000
                  -0.004110386631    -0.001612612835     0.000000000000
                  -0.000403748804     0.002889235801    -0.003063294403
                  -0.000403748804     0.002889235801     0.003063294403
                   0.001893737577    -0.002449099587    -0.002868342407
                   0.001893737577    -0.002449099587     0.002868342407],
    "formaldehyde" => [ 0.000000000000     0.000000000000     0.038875530921
                        0.000000000000     0.000000000000    -0.032184078783
                       -0.000000000000     0.006228089378    -0.003345726070
                        0.000000000000    -0.006228089378    -0.003345726070], 
    "glycine" => [-0.010904846656     0.017275971537     0.000000000000
                  -0.026799487548    -0.008060718767     0.000000000000
                   0.037341615244    -0.008058958568     0.000000000000
                   0.042776865478     0.002700985070     0.000000000000
                  -0.006721087456     0.016051709587     0.000000000000
                  -0.010497051729    -0.021349843046     0.000000000000
                  -0.008374522485     0.001083021282     0.007877987928
                  -0.008374522485     0.001083021282    -0.007877987928
                  -0.004223481179    -0.000362594189    -0.009527544824
                  -0.004223481179    -0.000362594189     0.009527544824],
    "methane" => [-0.000000000000    -0.000000000000     0.000000000000
                   0.001502203384     0.001502203384     0.001502203384
                  -0.001502203384    -0.001502203384     0.001502203384
                  -0.001502203384     0.001502203384    -0.001502203384
                   0.001502203384    -0.001502203384    -0.001502203384], 
    "phosphaethene" => [-0.004083749858    -0.010319837126     0.000000000000
                         0.000107570281     0.000056764675     0.000000000000
                        -0.006683153499     0.004032387469     0.000000000000
                         0.007800837937     0.002942709443     0.000000000000
                         0.002858495139     0.003287975562     0.000000000000], 
    "acetic_acid" => [-0.001515921968    -0.002091844521     0.000000000000
                      -0.028144294562     0.026907802176     0.000000000000
                       0.004462025923    -0.010008659808     0.000000000000
                      -0.000255682949     0.004993941401     0.000000000000
                      -0.000326139336    -0.002764392697     0.003360578328
                      -0.000326139336    -0.002764392697    -0.003360578328
                       0.009572222174    -0.029043547883     0.000000000000
                       0.016533930057     0.014771094030     0.000000000000], 
    "nitrogen" => [0.000000000000     0.000000000000     0.084169980572
                   0.000000000000     0.000000000000    -0.084169980572]
)

@testset "RHF" begin

    # Test argument error
    @test_throws Fermi.Options.FermiException Fermi.HartreeFock.RHF(1)

    # Test core guess 
    mol = Molecule()
    @set scf_guess core
    wfn = @energy mol => rhf
    @test isapprox(wfn.energy, -74.96500289468489, rtol=tol)

    # Test dense array
    Iu = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())
    wfn = @energy Iu => rhf
    @test isapprox(wfn.energy, -74.96500289468489, rtol=tol)

    # Test using initial wave function guess
    @set scf_max_iter 5
    wfn2 = @energy wfn => rhf
    @test isapprox(wfn.energy, -74.96500289468489, rtol=tol)

    # Test invalid number of electrons
    @set charge 1
    @set multiplicity 2
    mol = Molecule()
    @test_throws Fermi.Options.FermiException @energy mol => rhf

    @reset
    @set printstyle none
    @testset "Conventional" begin
        Fermi.Options.set("df", false)

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            Iu = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.SparseERI())

            wf = @energy Iu => rhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Density Fitted" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvqz-jkfit")

        for i = eachindex(molecules)
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("rhf_alg", 1)
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", basis[i])

            wf = @energy rhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4

            if i in [3,5,7]
                Fermi.Options.set("rhf_alg", 2)
                wf = @energy rhf
                @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
            end
        end
    end

    @testset "Gradients" begin
        @reset
        @set printstyle none

        # Test argument error         
        Fermi.Options.set("deriv_type", "Test error")
        @test_throws Fermi.Options.FermiException Fermi.HartreeFock.RHFgrad()

        Fermi.Options.set("deriv_type", "analytic")

        # Default gradient call
        for i = 1:10
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*molecules[i]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring",mol)
            Fermi.Options.set("basis", basis[i])

            g = @gradient rhf
            Δg = g - ref_grads[molecules[i]]
            @test √(sum(Δg.^2)/length(g)) < tol # Gradient from Psi4
        end

        # Chonky integrals
        @testset "Chonky" begin
            # Read molecule
            path = joinpath(@__DIR__, "xyz/water.xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring",mol)
            Fermi.Options.set("basis", "cc-pvtz")

            g = @gradient rhf <= Fermi.Integrals.Chonky()
            Δg = g - ref_grads["water"]
            @test √(sum(Δg.^2)/length(g)) < tol # Gradient from Psi4
        end
    end
end
