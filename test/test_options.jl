
# Simple test for macros

@testset "Options" begin
    # Basic check
    @set scf_e_conv 10^-10
    x = @get scf_e_conv
    @test x == 10^-10

    @set basis cc-pvqz
    x = @get basis
    @test x == "cc-pvqz"

    @set basis "6-31g*"
    x = @get basis
    @test x == "6-31g*"

    # Check fetching default
    y = @get cc_max_iter
    @test y == Fermi.Options.Default["cc_max_iter"]

    # Check error due to invalid key
    @test_throws Fermi.InvalidFermiOption @get invalid_options
    @test_throws Fermi.InvalidFermiOption @set invalid_option 10

    # Check error due to invalid data type
    @test_throws Fermi.InvalidFermiOption @set molstring 13

    # Check attempt to convert data type
    @set scf_max_iter 10.0
    x = @get scf_max_iter 
    @test x == 10
    @test_throws Fermi.InvalidFermiOption @set cc_max_iter 13.567

    # Check lowercase call (i.e. all options are saved in lowercase)
    @set {
        precision DOUBLE
    }
    x = @get precision
    @test x == "double"


    # Check if reset works
    @reset
    @test Fermi.Options.Current == Dict()
    @set {
        scf_e_conv 1/10
        cc_max_iter 1
        diis false
    }


    @reset diis
    @test @get diis 
    @reset cc_max_iter scf_e_conv
    x = @get scf_e_conv
    y = @get cc_max_iter
    @test x == Fermi.Options.Default["scf_e_conv"] && y == Fermi.Options.Default["cc_max_iter"]
    @test_throws Fermi.InvalidFermiOption @reset invalid_keyword 

    @reset
    @set {
        basis cc-pvdz
    }

    # @get (no input)
    x = @capture_out @get
    @test x == """
            ┌─────────┬───────────────┐
            │ Keyword │ Current Value │
            ├─────────┼───────────────┤
            │   basis │       cc-pvdz │
            └─────────┴───────────────┘\n"""

    # Check lookup
    @reset
    @set {
        scf_max_rms 1e-9
        scf_max_iter 50
    }

    x = @capture_out @lookup scf_max
    @test x == """
            ┌──────────────┬───────────────┐
            │      Keyword │ Current Value │
            ├──────────────┼───────────────┤
            │  scf_max_rms │   1.00000e-09 │
            │ scf_max_iter │            50 │
            └──────────────┴───────────────┘\n"""

end