using Molecules

@reset
@set printstyle none

Econv = [-75.4196061687054,  -39.5638172309568,  -38.9378598408943]
Edf =   [-75.41960080317435, -39.56380446089254, -38.93785583198393]
@testset "UHF" begin
    @testset "Conventional Chonky" begin
        Fermi.Options.set("df", false)
        for i = [1]
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            I = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())

            wf = @energy I => uhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end
    
    @testset "Conventional SparseERI" begin
        Fermi.Options.set("df", false)
        for i = [1]
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            I = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.SparseERI())

            wf = @energy I => uhf
            @test isapprox(wf.energy, Econv[i], rtol=tol) # Energy from Psi4
        end
    end

    @testset "Density Fitting" begin
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvtz-jkfit")
        for i in [1]
            # Read molecule
            path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
            mol = open(f->read(f,String), path)

            # Define options
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", uhf_boys[i][2])
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", uhf_boys[i][3])
            Fermi.Options.set("multiplicity", uhf_boys[i][4])
            Fermi.Options.set("uhf_alg", 1)

            wf = @energy uhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4

            # Direct (integral-direct DF) should agree with the in-core DF energy above
            Fermi.Options.set("uhf_alg", 2)
            wf = @energy uhf
            @test isapprox(wf.energy, Edf[i], rtol=tol) # Energy from Psi4
            Fermi.Options.set("uhf_alg", 1)
        end
    end

    @reset
    @set printstyle none
    @testset "Gradients" begin

        # Test argument error
        Fermi.Options.set("deriv_type", "Test error")
        @test_throws Fermi.Options.FermiException Fermi.HartreeFock.UHFgrad()

        Fermi.Options.set("deriv_type", "analytic")

        # Closed-shell reduction: UHF forced on a closed-shell singlet must
        # reduce to Dα=Dβ and reproduce the already-validated RHF gradient
        # to near machine precision -- see UHFgrad.jl's header comment for
        # the derivation this checks (X^A_UHF reduces exactly to RHF's X^A
        # in this limit). A much stronger, external-reference-independent
        # check than finite difference or Psi4 alone.
        @testset "Closed-shell reduction to RHF" begin
            Fermi.Options.set("molstring", """
            O        0.000000000000     0.000000000000     0.000000000000
            H        0.000000000000     0.000000000000     0.950000000000
            H        0.895669800000     0.000000000000    -0.316663000000
            """)
            Fermi.Options.set("basis", "cc-pvdz")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 1)

            rhf_wfn = @energy rhf
            rhf_g = Fermi.HartreeFock.RHFgrad(rhf_wfn)

            uhf_wfn = Fermi.HartreeFock.UHF()
            uhf_g = Fermi.HartreeFock.UHFgrad(uhf_wfn)

            @test isapprox(uhf_wfn.energy, rhf_wfn.energy, atol=1e-10)
            @test maximum(abs.(uhf_g .- rhf_g)) < 1e-9
        end

        # Open-shell references from a fresh Psi4 UHF gradient (scf_type
        # pk, matching Fermi's exact-ERI dispatch), same uhf_boys systems
        # the energy tests above already validate against Psi4 -- one
        # doublet (oh), one doublet radical (methyl), one triplet
        # (methylene), covering both Nα>Nβ>0 cases used elsewhere in this
        # file plus a genuinely different spin multiplicity.
        ref_uhf_grads = Dict(
            "oh" => [0.0 0.0 0.00012086175; 0.0 0.0 -0.00012086175],
            "methyl" => [-9.3854145e-5 -4.3407006e-5 -1.119603e-6;
                          1.3629195e-5  2.0166761e-5 -1.842176e-6;
                          2.26684e-6    9.317891e-6   7.38686e-7;
                          7.795811e-5   1.3922353e-5  2.223093e-6],
            "methylene" => [-0.0  0.0           9.737248e-6;
                             -0.0  7.1956668e-5 -4.868624e-6;
                              0.0 -7.1956668e-5 -4.868624e-6],
        )
        @testset "Psi4 reference (exact-ERI)" begin
            for i = eachindex(uhf_boys)
                path = joinpath(@__DIR__, "xyz/"*uhf_boys[i][1]*".xyz")
                mol = open(f->read(f,String), path)

                Fermi.Options.set("df", false)
                Fermi.Options.set("molstring", mol)
                Fermi.Options.set("basis", uhf_boys[i][2])
                Fermi.Options.set("reference", "uhf")
                Fermi.Options.set("charge", uhf_boys[i][3])
                Fermi.Options.set("multiplicity", uhf_boys[i][4])

                g = @gradient uhf
                ref = ref_uhf_grads[uhf_boys[i][1]]
                rms = sqrt(sum((g .- ref).^2)/length(g))
                @test rms < 1e-7
            end
        end

        # Nβ=0 edge case (untested by any of the uhf_boys systems above,
        # all of which have Nβ>=3): a bare doublet H atom. Regression test
        # for the odadamping 0/0 (NaN) bug this uncovered -- see
        # UHFHelper.jl's odadamping, dD is identically zero for the whole
        # SCF when a spin channel is empty, so c=tr((F-Fs)*dD)=0 exactly,
        # and the un-guarded -s/(2*c) branch condition was a literal 0/0.
        # Gradient must be exactly zero (single atom, translational
        # invariance).
        @testset "Nβ=0 edge case" begin
            Fermi.Options.set("df", false)
            Fermi.Options.set("molstring", "H 0.0 0.0 0.0")
            Fermi.Options.set("basis", "cc-pvdz")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 2)

            wfn = Fermi.HartreeFock.UHF()
            @test isfinite(wfn.energy)
            g = Fermi.HartreeFock.UHFgrad(wfn)
            @test all(iszero, g)
        end
    end

    @reset
    @set printstyle none
    @testset "Hessian" begin

        # Closed-shell reduction: same rationale as the Gradients testset's
        # version -- UHF forced on a closed-shell singlet must reduce to
        # Dα=Dβ and reproduce the already-validated RHF Hessian to near
        # machine precision. Also the check that caught a real bug during
        # development: UCPHF's Fskel used RHF's CPHF.jl's `2*Jq` coefficient
        # verbatim, but UHF's Fα=H+Jtot-Kα already has Jtot built directly
        # from Dtot (no extra factor belongs there) -- this test's diff went
        # from ~2.4 to ~2e-11 once that was fixed.
        @testset "Closed-shell reduction to RHF" begin
            Fermi.Options.set("molstring", """
            O        0.000000000000     0.000000000000     0.000000000000
            H        0.000000000000     0.000000000000     0.950000000000
            H        0.895669800000     0.000000000000    -0.316663000000
            """)
            Fermi.Options.set("basis", "sto-3g")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 1)

            rhf_wfn = @energy rhf
            rhf_H = Fermi.HartreeFock.RHFhess(rhf_wfn)

            uhf_wfn = Fermi.HartreeFock.UHF()
            uhf_H = Fermi.HartreeFock.UHFhess(uhf_wfn)

            @test uhf_H ≈ uhf_H' atol=1e-9
            @test maximum(abs.(uhf_H .- rhf_H)) < 1e-8
        end

        # Independent, open-shell check: finite difference of Fermi's own
        # analytic UHF gradient (already validated above), same pattern
        # test_hessian.jl's RHF version uses.
        @testset "Finite difference (open-shell)" begin
            Fermi.Options.set("molstring", """
            O   0.000000000000   0.000000000000   0.000000000000
            H   0.000000000000   0.000000000000   0.970000000000
            """)
            basis = "sto-3g"
            Fermi.Options.set("basis", basis)
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 2)

            wf = Fermi.HartreeFock.UHF()
            H = Fermi.HartreeFock.UHFhess(wf)
            @test H ≈ H' atol=1e-10

            function displaced_molstring(atoms, iA, k, delta)
                newatoms = deepcopy(atoms)
                xyz = collect(newatoms[iA].xyz)
                xyz[k] += delta
                newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
                return join(["$(Int(a.Z))   $(a.xyz[1])   $(a.xyz[2])   $(a.xyz[3])" for a in newatoms], "\n")
            end

            function uhf_grad(molstring)
                Fermi.Options.set("molstring", molstring)
                Fermi.Options.set("basis", basis)
                Fermi.Options.set("charge", 0)
                Fermi.Options.set("multiplicity", 2)
                wf_ = Fermi.HartreeFock.UHF()
                return Fermi.HartreeFock.UHFgrad(wf_)
            end

            atoms = wf.molecule.atoms
            natm = length(atoms)
            B2A = Molecules.bohr_to_angstrom
            h = 5e-4
            H_FD = zeros(3*natm, 3*natm)
            for iB in 1:natm, qB in 1:3
                gp = uhf_grad(displaced_molstring(atoms, iB, qB, h))
                gm = uhf_grad(displaced_molstring(atoms, iB, qB, -h))
                H_FD[:, 3*(iB-1)+qB] .= vec(((gp .- gm) ./ (2h/B2A))')
            end
            @test H ≈ H_FD atol=1e-5
        end

        # Fresh Psi4 UHF Hessian reference (scf_type pk, matching Fermi's
        # exact-ERI dispatch), oh/cc-pvtz -- same system the gradient's
        # Psi4-reference testset uses.
        @testset "Psi4 reference (exact-ERI)" begin
            path = joinpath(@__DIR__, "xyz/oh.xyz")
            mol = open(f->read(f,String), path)
            Fermi.Options.set("df", false)
            Fermi.Options.set("molstring", mol)
            Fermi.Options.set("basis", "cc-pvtz")
            Fermi.Options.set("reference", "uhf")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 2)

            wfn = Fermi.HartreeFock.UHF()
            H = Fermi.HartreeFock.UHFhess(wfn)

            H_psi4 = [
             -0.00006725075105  0.0                 0.0                 0.00006725075183   0.0                 0.0
              0.0               -0.00006724582223    0.0                 0.0                 0.00006724582304    0.0
              0.0                0.0                 0.58571680342079    0.0                 0.0                -0.58571680342115
              0.00006725075183   0.0                 0.0                -0.00006725075189   0.0                 0.0
              0.0                0.00006724582304    0.0                 0.0                -0.00006724582306    0.0
              0.0                0.0                -0.58571680342115    0.0                 0.0                 0.58571680342118
            ]
            @test maximum(abs.(H .- H_psi4)) < 1e-7
        end

        # Nβ=0 edge case: bare doublet H atom. Hessian must be exactly zero
        # (single atom, translational invariance) and finite (no NaN from
        # the odadamping bug the gradient testset already regression-tests,
        # or from any UCPHF-side empty-spin-channel edge case).
        @testset "Nβ=0 edge case" begin
            Fermi.Options.set("df", false)
            Fermi.Options.set("molstring", "H 0.0 0.0 0.0")
            Fermi.Options.set("basis", "cc-pvdz")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 2)

            wfn = Fermi.HartreeFock.UHF()
            H = Fermi.HartreeFock.UHFhess(wfn)
            @test all(x -> isfinite(x) && abs(x) < 1e-8, H)
        end

        @testset "@hessian / @frequencies macro dispatch" begin
            Fermi.Options.set("molstring", """
            O        0.000000000000     0.000000000000     0.000000000000
            H        0.000000000000     0.000000000000     0.950000000000
            H        0.895669800000     0.000000000000    -0.316663000000
            """)
            Fermi.Options.set("basis", "sto-3g")
            Fermi.Options.set("charge", 0)
            Fermi.Options.set("multiplicity", 1)

            wfn = Fermi.HartreeFock.UHF()
            H_direct = Fermi.HartreeFock.UHFhess(wfn)
            H_macro = @hessian wfn => uhf
            @test H_macro ≈ H_direct atol=1e-12

            freqs_ref, modes_ref = vibrational_analysis(H_direct, wfn.molecule)
            freqs, modes = @frequencies wfn => uhf
            @test freqs ≈ freqs_ref atol=1e-8
            @test modes ≈ modes_ref atol=1e-8

            # @frequency is an alias for @frequencies -- same result
            freqs2, modes2 = @frequency wfn => uhf
            @test freqs2 ≈ freqs_ref atol=1e-8
            @test modes2 ≈ modes_ref atol=1e-8
        end
    end
end