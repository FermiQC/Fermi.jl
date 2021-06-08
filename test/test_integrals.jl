@testset "Integrals" begin
    @reset

    @set {
        jkfit auto
        rifit auto
        basis cc-pvtz
    }

    jkfit = Fermi.Integrals.JKFIT()
    rifit = Fermi.Integrals.RIFIT()

    # Test auto-complete
    @test jkfit.basisset.basis_name == "cc-pvtz-jkfit"
    @test rifit.basisset.basis_name == "cc-pvtz-rifit"

    # Test displays
    x = @capture_out display(jkfit)
    @test occursin(r"JKFIT:\s+?cc-pvtz-jkfit",x)
    x = @capture_out display(rifit)
    @test occursin(r"RIFIT:\s+?cc-pvtz-rifit",x)

    ints = Fermi.Integrals.IntegralHelper()

    # Try invalid precision
    @set precision invalid
    @test_throws Fermi.Options.FermiException Fermi.Integrals.IntegralHelper()

    @set precision single
    ints32 = Fermi.Integrals.IntegralHelper()

    # Test if single and double are approx equal
    @test isapprox(ints["S"], ints32["S"])
    @test isapprox(ints["T"], ints32["T"])
    @test isapprox(ints["V"], ints32["V"])
    @test isapprox(ints["ERI"], ints32["ERI"])

    # Test string_repr
    x = @capture_out display(ints)
    @test begin
        occursin(r"Fermi\s+?IntegralHelper",x) &&
        occursin(r"Data\s+?Type:\s+?Float64",x) &&
        occursin(r"Basis:\s+?cc-pvtz",x) &&
        occursin(r"ERI:\s+?Chonky",x) &&
        occursin(r"Orbitals:\s+?AtomicOrbitals",x)
    end

    # Test cache
    @test begin
        "S" in keys(ints.cache) &&
        "T" in keys(ints.cache) &&
        "V" in keys(ints.cache) &&
        "ERI" in keys(ints.cache)
    end

    # Test two basis set input
    bf = ints.orbitals.basisset
    S = Fermi.Integrals.ao_1e(bf, bf, "overlap")
    T = Fermi.Integrals.ao_1e(bf, bf, "kinetic")
    V = Fermi.Integrals.ao_1e(bf, bf, "nuclear")

    @test S ≈ ints["S"] 
    @test T ≈ ints["T"] 
    @test V ≈ ints["V"]

    # Test delete one key
    delete!(ints, "V")
    @test begin
        "S" in keys(ints.cache) &&
        "T" in keys(ints.cache) &&
      !("V" in keys(ints.cache)) &&
        "ERI" in keys(ints.cache)
    end
    # Test deleting everything
    delete!(ints)
    @test isempty(ints.cache)

    # Test invalid key
    @test_throws Fermi.Options.FermiException ints["invalid"]

end