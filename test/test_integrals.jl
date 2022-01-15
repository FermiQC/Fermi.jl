@testset "Misc Integrals" begin
    @reset
    aoints = Fermi.Integrals.IntegralHelper()
    @test_throws Fermi.Options.FermiException Fermi.Integrals.compute!(aoints, "invalid")

    @set {
        basis = cc-pvtz
        jkfit auto
        rifit auto
    }

    jkfit = Fermi.Integrals.JKFIT()
    @test Fermi.string_repr(jkfit) == "JKFIT: cc-pvtz-jkfit"

    rifit = Fermi.Integrals.RIFIT()
    @test Fermi.string_repr(rifit) == "RIFIT: cc-pvtz-rifit"

    @reset
    @set precision "xablau"
    @test_throws Fermi.Options.FermiException Fermi.Integrals.IntegralHelper()

    @set precision "single"
    ints = Fermi.Integrals.IntegralHelper()
    S = ints["S"]
    @test S[1,1] |> typeof === Float32

    ints["T"]
    ints["V"] = rand(5,5)
    @test begin
        "S" in keys(ints.cache) && 
        "T" in keys(ints.cache) &&
        "V" in keys(ints.cache)
    end

    Fermi.Integrals.delete!(ints, "V")
    @test begin
        "S" in keys(ints.cache) && 
        "T" in keys(ints.cache) &&
        !("V" in keys(ints.cache))
    end
    Fermi.Integrals.delete!(ints)

    @test begin
        !("S" in keys(ints.cache)) && 
        !("T" in keys(ints.cache)) &&
        !("V" in keys(ints.cache))
    end

    @reset
    mol = Molecule(molstring="He 0.0 0.0 0.0")
    bset = Fermi.Integrals.BasisSet("sto-3g", mol)
    ints = Fermi.Integrals.IntegralHelper(bset)
    ints["S"]
    @test Fermi.string_repr(ints) == " ⇒ Fermi IntegralHelper\n ⋅ Data Type:                 Float64\n ⋅ Basis:                     sto-3g\n ⋅ ERI:                       SparseERI\n ⋅ Orbitals:                  AtomicOrbitals\n ⋅ Stored Integrals:          S "

end

@testset "1-e AO -> MO" begin
    @reset
    @set {
        printstyle none
        df true
        basis cc-pvdz
        jkfit cc-pvtz-jkfit
    }
    aoints = Fermi.Integrals.IntegralHelper()
    wfn = @energy aoints => rhf
    Cμi = wfn.orbitals.C
    moints = Fermi.Integrals.IntegralHelper(orbitals=wfn.orbitals)

    Fermi.Integrals.compute!(moints, aoints, "S")
    @test moints.cache["S"] ≈ Cμi' * aoints["S"] * Cμi

    Fermi.Integrals.compute!(moints, aoints, "T")
    @test moints.cache["T"] ≈ Cμi' * aoints["T"] * Cμi

    Fermi.Integrals.compute!(moints, aoints, "V")
    @test moints.cache["V"] ≈ Cμi' * aoints["V"] * Cμi
end

@testset "DF AO -> MO" begin

    @reset
    @set {
        printstyle none
        df true
        basis cc-pvdz
        jkfit cc-pvtz-jkfit
        drop_occ 2
    }

    @molecule {
        C       -1.1170269915      4.4925791088      0.0000000582                 
        C        0.3969669792      4.3675789789     -0.0000001928                 
        H       -1.5729510532      3.5466015111      0.3608481999                 
        H       -1.4776736551      4.7005983883     -1.0293753400                 
        H       -1.4262457918      5.3234861671      0.6685275007                 
        H        0.7061857730      3.5366719134     -0.6685276289                 
        H        0.7576136432      4.1595597057      1.0293752066                 
        H        0.8528910462      5.3135565710     -0.3608483431                 

    }

    wfn = @energy rhf
    Cμi = wfn.orbitals.C
    Co = Cμi[:, 3:9]
    Cv = Cμi[:,10:end]

    moints = Fermi.Integrals.IntegralHelper(orbitals=wfn.orbitals, eri_type=Fermi.Integrals.RIFIT())
    aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.RIFIT())

    # contract B(prs) = C(μr) * B(pμν) * C(νs) from scratch
    function test_contract(C1, C2, B)
        p,μ,ν = size(B)
        r = size(C1, 2)
        s = size(C2, 2)

        Bμpν = permutedims(B, (2,1,3))
        Bμpν = reshape(Bμpν, (μ, p*ν))
        Brpν = C1' * Bμpν 
        Brpν = reshape(Brpν, (r*p, ν))
        Brps = Brpν * C2
        Brps = reshape(Brps, (r,p,s))
        return permutedims(Brps, (2,1,3))
    end

    # Contract (ijkl) = B(pij) * B(pkl)
    function test_contract2(B1, B2)
        p,i,j = size(B1)
        _,k,l = size(B2)
        B1 = reshape(B1, p,i*j)
        B2 = reshape(B2, p,k*l)
        out = B1' * B2
        return reshape(out, (i,j,k,l))
    end

    # Full Bpq
    @test moints["ERI"] ≈ test_contract(Cμi, Cμi, aoints["ERI"])

    # Boo
    @test moints["BOO"] ≈ test_contract(Co, Co, aoints["ERI"])

    # Bov
    @test moints["BOV"] ≈ test_contract(Co, Cv, aoints["ERI"])

    # Bvv
    @test moints["BVV"] ≈ test_contract(Cv, Cv, aoints["ERI"])

    Boo = moints["BOO"]
    Bov = moints["BOV"]
    Bvv = moints["BVV"]

    oooo = test_contract2(Boo, Boo)
    ooov = test_contract2(Boo, Bov)
    oovv = test_contract2(Boo, Bvv)
    ovov = test_contract2(Bov, Bov)
    ovvv = test_contract2(Bov, Bvv)
    vvvv = test_contract2(Bvv, Bvv)

    @test moints["OOOO"] ≈ oooo
    @test moints["OOOV"] ≈ ooov
    @test moints["OOVV"] ≈ oovv
    @test moints["OVOV"] ≈ ovov
    @test moints["OVVV"] ≈ ovvv
    @test moints["VVVV"] ≈ vvvv

    # Test again, cleaning the cache such that Bpq arrays need to be compute on the fly
    delete!(moints)

    @test moints["OOOO"] ≈ oooo
    @test moints["OOOV"] ≈ ooov
    @test moints["OOVV"] ≈ oovv
    @test moints["OVOV"] ≈ ovov
    @test moints["OVVV"] ≈ ovvv
    @test moints["VVVV"] ≈ vvvv
end


@testset "Chonky AO -> MO" begin

    @reset
    @set {
        printstyle none
        basis cc-pvdz
        drop_occ 1
    }

    @molecule {
        C        0.3340594014      4.3347983152      0.0216449677                 
        H        1.4398420378      4.2538788474     -0.0166555421                 
        H       -0.0088087137      5.0751637968     -0.7300667682                 
        H        0.0222217713      4.6642047575      1.0340770936                 
        H       -0.1170174852      3.3459458686     -0.2007749095                 
    }

    wfn = @energy rhf
    C  = wfn.orbitals.C
    Co = C[:, 2:5]
    Cv = C[:,6:end]

    moints = Fermi.Integrals.IntegralHelper(orbitals=wfn.orbitals, eri_type=Fermi.Integrals.Chonky())
    aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())

    # contract pqrs = C1(μp) * C2(νq) * C3(ρr) * C4(σs) * μνρσ
    function test_contract(C1, C2, C3, C4, μνρσ)

        μ,ν,ρ,σ = size(μνρσ)
        _,s = size(C4)
        _,r = size(C3)
        _,q = size(C2)
        _,p = size(C1)

        # μνρ,σ * σs
        μνρs = reshape(reshape(μνρσ, (μ*ν*ρ, σ)) * C4, (μ,ν,ρ,s))

        # μνs,ρ * ρr 
        μνsr = reshape(reshape(permutedims(μνρs, (1,2,4,3)), (μ*ν*s, ρ)) * C3, (μ,ν,s,r))

        # μrs,ν * νq
        μrsq = reshape(reshape(permutedims(μνsr, (1,4,3,2)), (μ*r*s, ν)) * C2, (μ,r,s,q))

        # rsq,μ * μp
        pqrs = reshape(C1' * reshape(permutedims(μrsq, (1,4,2,3)), (μ, q*r*s)), (p,q,r,s))

        return pqrs
    end

    # Full ERI
    @test moints["ERI"] ≈ test_contract(C, C, C, C, aoints["ERI"])
    delete!(moints)

    # oooo
    @test moints["OOOO"] ≈ test_contract(Co, Co, Co, Co, aoints["ERI"])

    # ooov
    @test moints["OOOV"] ≈ test_contract(Co, Co, Co, Cv, aoints["ERI"])

    # oovv
    @test moints["OOVV"] ≈ test_contract(Co, Co, Cv, Cv, aoints["ERI"])

    # ovov
    @test moints["OVOV"] ≈ test_contract(Co, Cv, Co, Cv, aoints["ERI"])

    # ovvv
    @test moints["OVVV"] ≈ test_contract(Co, Cv, Cv, Cv, aoints["ERI"])

    # vvvv
    @test moints["VVVV"] ≈ test_contract(Cv, Cv, Cv, Cv, aoints["ERI"])
end

@testset "Sparse AO -> MO" begin

    @reset
    @set {
        printstyle none
        basis cc-pvdz
        drop_occ 1
        drop_vir 3
    }

    @molecule {
        N       -1.0811450850      4.5346506907     -0.6162384680                 
        H       -0.0380366676      4.5141429777     -0.6171266138                 
        H       -1.4191680136      3.9897206949      0.2067368203                 
        H       -1.4191724649      4.0628676478     -1.4832192333                 
    }

    wfn = @energy rhf

    mo1 = Fermi.Integrals.IntegralHelper(orbitals=wfn.orbitals, eri_type=Fermi.Integrals.Chonky())
    mo2 = Fermi.Integrals.IntegralHelper(orbitals=wfn.orbitals, eri_type=Fermi.Integrals.Chonky())
    aosp = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.SparseERI())
    aock = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())

    # oooo
    Fermi.Integrals.compute_OOOO!(mo1, aosp)
    Fermi.Integrals.compute_OOOO!(mo2, aock)
    @test mo1["OOOO"] ≈ mo2["OOOO"]

    # ooov
    Fermi.Integrals.compute_OOOV!(mo1, aosp)
    Fermi.Integrals.compute_OOOV!(mo2, aock)
    @test mo1["OOOV"] ≈ mo2["OOOV"]

    # oovv
    Fermi.Integrals.compute_OOVV!(mo1, aosp)
    Fermi.Integrals.compute_OOVV!(mo2, aock)
    @test mo1["OOVV"] ≈ mo2["OOVV"]

    # ovov
    Fermi.Integrals.compute_OVOV!(mo1, aosp)
    Fermi.Integrals.compute_OVOV!(mo2, aock)
    @test mo1["OVOV"] ≈ mo2["OVOV"]

    # ovvv
    Fermi.Integrals.compute_OVVV!(mo1, aosp)
    Fermi.Integrals.compute_OVVV!(mo2, aock)
    @test mo1["OVVV"] ≈ mo2["OVVV"]

    # vvvv
    Fermi.Integrals.compute_VVVV!(mo1, aosp)
    Fermi.Integrals.compute_VVVV!(mo2, aock)
    @test mo1["VVVV"] ≈ mo2["VVVV"]
end

