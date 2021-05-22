@testset "Basis Set" begin
    bfH = Fermi.GaussianBasis.BasisFunction(0, 
    [0.9817067282501578, 0.9494640079478909, 0.2959064596969707], 
    [3.425250914, 0.6239137298, 0.168855404])

    # Read a simple basis function
    x = Fermi.GaussianBasis.read_basisset("sto-3g", "H")[1] 
    @test begin
        x.l == bfH.l &&
        x.coef == bfH.coef &&
        x.exp == bfH.exp
    end

    # Test error: Invalid Basis set
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.read_basisset("invalid_basis", "H")

    # Test error: Invalid Atom
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.read_basisset("sto-3g", "U")

    # Test some basis set non sense
    x = """W    3   1.00
    0.3425250914D+01       0.1543289673D+0
    0.6239137298D+00       0.5353281423D+00
    0.1688554040D+00       0.4446345422D+00"""
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.basis_from_string(x)
    x = """S    3   1.00
    0.3425250914D+01       0.1543289673D+0
    0.6239137298D+00       Ops 
    0.1688554040D+00       0.4446345422D+00"""
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.basis_from_string(x)

    x =  """SP D   3   1.00
          0.5033151319D+01      -0.9996722919D-01       0.1559162750D+00
          0.1169596125D+01       0.3995128261D+00       0.6076837186D+00
          0.3803889600D+00       0.7001154689D+00       0.3919573931D+00"""
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.two_basis_from_string(x)

    x =  """Some   3   1.00
          0.5033151319D+01      -0.9996722919D-01       0.1559162750D+00
          0.1169596125D+01       0.3995128261D+00       0.6076837186D+00
          0.3803889600D+00       0.7001154689D+00       0.3919573931D+00"""
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.two_basis_from_string(x)

    x =  """SP   3   1.00
          0.5033151319D+01      -0.9996722919D-01       0.1559162750D+00
          0000000000000000       0.3995128261D+00       0.6076837186D+00
          0.3803889600D+00       0.7001154689D+00       0.3919573931D+00"""
    @test_throws Fermi.Options.FermiException Fermi.GaussianBasis.two_basis_from_string(x)

    @set basis sto-3g

    bs = Fermi.GaussianBasis.BasisSet()
    x = @capture_out display(bs)
    @test begin
          occursin("sto-3g Basis Set", x) &&
          occursin(r"Number of shells:\s+?5", x) &&
          occursin(r"Number of basis:\s+?7", x) &&
          occursin("O: 1s 2s 1p", x) &&
          occursin("H: 1s", x)
    end
    x = @capture_out display(bs[1,1])
    @test begin
          occursin("S shell with 1 basis built from 3 primitive gaussians", x) &&
          occursin(r"χ₀₀\s+?=\s+?\d+?\.\d+?⋅Y₀₀⋅exp\(-\d+?\.\d+?⋅r²\)", x) &&
          occursin(r"[+-]\s+?\d+?\.\d+?⋅Y₀₀⋅exp\(-\d+?\.\d+?⋅r²\)", x)
    end
    x = @capture_out display(bs[1][3])
    @test begin
          occursin("P shell with 3 basis built from 3 primitive gaussians", x) &&
          occursin(r"χ.{2,3}\s+?=\s+?\d+?\.\d+?⋅Y.{2,3}⋅r¹⋅exp\(-\d+?\.\d+?⋅r²\)", x) &&
          occursin(r"[+-]\s+?\d+?\.\d+?⋅Y.{2,3}⋅r¹⋅exp\(-\d+?\.\d+?⋅r²\)", x)
    end
end