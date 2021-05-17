using Fermi.Options

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
    @test y == Options.Default["cc_max_iter"]

    # Check error due to invalid key
    @test_throws FermiException @get invalid_options
    @test_throws FermiException @set invalid_option 10

    # Check error due to invalid data type
    @test_throws FermiException @set molstring 13

    # Check attempt to convert data type
    @set scf_max_iter 10.0
    x = @get scf_max_iter 
    @test x == 10
    @test_throws FermiException @set cc_max_iter 13.567

    # Check lowercase call (i.e. all options are saved in lowercase)
    @set {
        precision DOUBLE
    }
    x = @get precision
    @test x == "double"


    # Check if reset works
    @reset
    @test Options.Current == Dict()
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
    @test x == Options.Default["scf_e_conv"] && y == Options.Default["cc_max_iter"]
    @test_throws FermiException @reset invalid_keyword 

    @reset
    x = @capture_out @get
    @test strip(x) == "No user defined keywords found."

    @set {
        basis cc-pvdz
        diis false
        molstring "H 1 2 3\nH 3 4 5"
    }

    # @get (no input)
    x = @capture_out @get
    @test begin
        occursin("Keyword", x)       &&
        occursin("Current Value", x) &&
        occursin("basis", x)         &&
        occursin("diis", x)          &&
        occursin("false", x)         &&
        occursin("molstring", x)     &&
        occursin(r"H\s+?1\.??0*?\s+?2\.??0*?\s+?3\.??0*?", x)
    end

    # Check lookup
    @reset
    x = @capture_out @lookup someinvalidkeyword
    @test strip(x) == "No keywords found containing: someinvalidkeyword"

    @set {
        scf_max_rms 1e-9
        scf_max_iter 50
        molstring "H 1 2 3\nH 3 4 5"
    }

    x = @capture_out @lookup scf_max
    @test begin
        occursin("Keyword", x)       &&
        occursin("Current Value", x) &&
        occursin("scf_max_rms", x)   &&
        occursin("1.00000e-09", x)   &&
        occursin("scf_max_iter", x)  &&
        occursin("50", x)  
    end
    x = @capture_out @lookup molstring
    @test begin
        occursin("molstring", x)     &&
        occursin(r"H\s+?1\.??0*?\s+?2\.??0*?\s+?3\.??0*?", x)
    end

    # Test @molecule
    molstring = """
    H 1.0 2.0 3.0
    H 4.0 5.0 6.0"""

    @molecule {
    H 1.0 2.0 3.0
    H 4.0 5.0 6.0
    }
    @test begin
        x = Options.get("molstring")
        occursin(r"H\s+?1\.0*?\s+?2\.0*?\s+?3\.0*?", x) &&
        occursin(r"H\s+?4\.0*?\s+?5\.0*?\s+?6\.0*?", x)
    end
    x = @capture_out showerror(IOContext(stdout, :short => true), Fermi.Options.FermiException("Testing the error!"))
    @test occursin("Testing the error!", x)
end