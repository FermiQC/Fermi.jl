using GaussianBasis
using Molecules
using TensorOperations

@reset
@set printstyle none

# Validates the two-electron (Coulomb+exchange) contribution to the RHF Hessian,
# ERI_hess_JK, against finite difference of the *same* 0.5*P*P*(ERI - 0.5*ERI_exch)
# contraction RHFgrad.jl already uses at first order -- i.e. this checks the
# "integral response" piece only (P held fixed), not the full CPHF-coupled
# Hessian, which doesn't exist yet.
@testset "Hessian" begin
    @testset "Two-electron JK contraction (finite difference)" begin
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")

        wf = @energy rhf
        bset = GaussianBasis.BasisSet(wf.orbitals.basis, wf.molecule.atoms)
        Co = wf.orbitals.C[:, 1:wf.ndocc]
        P = 2.0 .* Co * Co'

        function E2_coulomb(bs, P)
            eri = GaussianBasis.ERI_2e4c(bs)
            return sum(P[m,n]*P[r,s]*eri[m,n,r,s] for m in 1:bs.nbas, n in 1:bs.nbas, r in 1:bs.nbas, s in 1:bs.nbas)
        end

        function E2_exchange(bs, P)
            eri = GaussianBasis.ERI_2e4c(bs)
            return sum(P[m,n]*P[r,s]*eri[m,s,r,n] for m in 1:bs.nbas, n in 1:bs.nbas, r in 1:bs.nbas, s in 1:bs.nbas)
        end

        function FD2(f, bs, iA, iB, h=1e-4)
            out = zeros(3,3)
            B2A = Molecules.bohr_to_angstrom
            for kA in 1:3, kB in 1:3
                bAp, bAm = GaussianBasis.create_displacement(bs, iA, kA, h)
                bApBp, _ = GaussianBasis.create_displacement(bAp, iB, kB, h)
                _, bApBm = GaussianBasis.create_displacement(bAp, iB, kB, h)
                bAmBp, _ = GaussianBasis.create_displacement(bAm, iB, kB, h)
                _, bAmBm = GaussianBasis.create_displacement(bAm, iB, kB, h)
                out[kA,kB] = (f(bApBp) - f(bApBm) - f(bAmBp) + f(bAmBm)) / (4*(h/B2A)^2)
            end
            return out
        end

        natoms = length(wf.molecule.atoms)
        for iA in 1:natoms, iB in iA:natoms
            H2e_analytic = Fermi.HartreeFock.ERI_hess_JK(bset, P, iA, iB)
            H2e_FD = 0.5 .* FD2(bs_->E2_coulomb(bs_,P), bset, iA, iB) .-
                     0.25 .* FD2(bs_->E2_exchange(bs_,P), bset, iA, iB)
            @test H2e_analytic ≈ H2e_FD atol=1e-4
        end
    end

    @testset "CPHF response-Fock matrix-vector product" begin
        # Validates cphf_Amatvec (the CPHF coupling matrix's action A*U) against
        # a brute-force dense MO-ERI contraction of sum_bj[4(ai|bj)-(ab|ij)-(aj|ib)]U_bj,
        # a structurally independent computation (explicit AO->MO transform of the
        # full ERI tensor, not the build_fock!-reuse trick cphf_Amatvec relies on).
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", "sto-3g")

        wf = @energy rhf
        ndocc = wf.ndocc
        nbas = size(wf.orbitals.C, 1)
        nvir = nbas - ndocc
        C = wf.orbitals.C
        Co = C[:, 1:ndocc]
        Cv = C[:, ndocc+1:end]

        ints = Fermi.Integrals.IntegralHelper(molecule=wf.molecule, eri_type=Fermi.Integrals.Chonky())
        eri_ao = ints["ERI"]

        function mo_eri(eri_ao, Ca, Cb, Cc, Cd)
            @tensoropt eri_mo[a,b,c,d] := Ca[m,a]*Cb[n,b]*Cc[r,c]*Cd[s,d]*eri_ao[m,n,r,s]
            return eri_mo
        end

        eri_aibj = mo_eri(eri_ao, Cv, Co, Cv, Co)
        eri_abij = mo_eri(eri_ao, Cv, Cv, Co, Co)
        eri_ajib = mo_eri(eri_ao, Cv, Co, Co, Cv)

        # deterministic pseudo-random trial matrix (avoids a Random stdlib
        # test-dependency for what only needs to be "generic, not special")
        U = [sin(1.3*a + 2.7*i) for a in 1:nvir, i in 1:ndocc]

        AU_brute = zeros(nvir, ndocc)
        for a in 1:nvir, i in 1:ndocc
            s = 0.0
            for b in 1:nvir, j in 1:ndocc
                s += (4*eri_aibj[a,i,b,j] - eri_abij[a,b,i,j] - eri_ajib[a,j,i,b]) * U[b,j]
            end
            AU_brute[a,i] = s
        end

        AU_fast = Fermi.HartreeFock.cphf_Amatvec(U, Co, Cv, ints)
        @test AU_fast ≈ AU_brute atol=1e-10
    end

    @testset "CPHF orbital response (finite difference of the SCF density)" begin
        # The definitive CPHF check: does the analytic dD/dy built from the
        # CPHF-solved U actually match how the *real*, re-converged SCF
        # density changes under a geometry displacement? This is independent
        # of every other test here -- it reruns full SCF at displaced
        # geometries rather than reusing any of cphf_rhs/cphf_Amatvec's own
        # machinery.
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)
        basis = "sto-3g"

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", basis)
        wf = @energy rhf

        ndocc = wf.ndocc
        C = wf.orbitals.C
        Co = C[:, 1:ndocc]
        Cv = C[:, ndocc+1:end]
        bset = GaussianBasis.BasisSet(wf.orbitals.basis, wf.molecule.atoms)
        ints = Fermi.Integrals.IntegralHelper(molecule=wf.molecule, eri_type=Fermi.Integrals.Chonky())
        atoms = wf.molecule.atoms

        function displaced_molstring(atoms, iA, k, delta)
            newatoms = deepcopy(atoms)
            xyz = collect(newatoms[iA].xyz)
            xyz[k] += delta
            newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
            return join(["$(Int(a.Z))   $(a.xyz[1])   $(a.xyz[2])   $(a.xyz[3])" for a in newatoms], "\n")
        end

        function scf_D(molstring)
            Fermi.Options.set("molstring", molstring)
            Fermi.Options.set("basis", basis)
            wf_ = @energy rhf
            Co_ = wf_.orbitals.C[:, 1:wf_.ndocc]
            return Co_ * Co_'
        end

        B2A = Molecules.bohr_to_angstrom
        h = 5e-4

        for iA in 1:length(atoms)
            U = Fermi.HartreeFock.cphf_solve(wf, ints, iA)

            ∂S = zeros(size(Co, 1), size(Co, 1), 3)
            GaussianBasis.∇overlap!(∂S, bset, iA)

            for q in 1:3
                Dp = scf_D(displaced_molstring(atoms, iA, q, h))
                Dm = scf_D(displaced_molstring(atoms, iA, q, -h))
                dD_FD = (Dp .- Dm) ./ (2h / B2A)

                ΔD = Cv * U[:, :, q] * Co'
                ΔD .+= ΔD'
                S_oo = Co' * (@view ∂S[:, :, q]) * Co
                dD_analytic = ΔD .- Co * S_oo * Co'

                @test dD_analytic ≈ dD_FD atol=1e-5
            end
        end

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", basis)
    end

    @testset "DF Jq/Kq and CPHF orbital response" begin
        # Phase 3 of the DF-RHF analytic Hessian plan: Jq_DF/Kq_DF (the
        # matrix-valued, density-fitted analog of eri_grad_JK) and the CPHF
        # machinery built on top of them (cphf_rhs_df/cphf_solve_df).
        Fermi.Options.set("molstring", """
        O   0.000000000000   0.000000000000  -0.143225816552
        H   0.000000000000   1.638036840407   1.136548822547
        H   0.000000000000  -1.638036840407   1.136548822547
        """)
        Fermi.Options.set("unit", "bohr")
        basis = "sto-3g"
        Fermi.Options.set("basis", basis)
        Fermi.Options.set("df", true)

        aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.RIFIT())
        wf = Fermi.HartreeFock.RHF(aoints, Fermi.HartreeFock.get_rhf_alg())
        cache = Fermi.HartreeFock.build_df_hess_cache(wf, aoints)
        atoms = wf.molecule.atoms
        natm = length(atoms)

        @testset "Kq_DF symmetry (free, FD-independent)" begin
            for iA in 1:natm
                _, Kq = Fermi.HartreeFock.eri_grad_JK_df(cache, iA)
                for q in 1:3
                    @test Kq[:, :, q] ≈ Kq[:, :, q]' atol=1e-8
                end
            end
        end

        @testset "Jq_DF/Kq_DF vs dense brute force" begin
            # Build the dense (mn|rs)_DF tensor from Pmn/Jinv directly and
            # finite-difference IT (a completely independent computation
            # path from Jq_DF/Kq_DF's own analytic formula), mirroring how
            # eri_grad_JK's own original validation worked one derivative
            # order down.
            bset, auxbset = cache.bset, cache.auxbset
            D = cache.D
            h = 1e-4
            B2A = Molecules.bohr_to_angstrom

            function displaced_bsets(iA, k, delta)
                newatoms = deepcopy(atoms)
                xyz = collect(newatoms[iA].xyz)
                xyz[k] += delta
                newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
                return BasisSet(basis, newatoms), BasisSet(aoints.eri_type.basisset.name, newatoms)
            end

            for iA in 1:natm
                Jq, Kq = Fermi.HartreeFock.eri_grad_JK_df(cache, iA)
                for k in 1:3
                    bsp, auxp = displaced_bsets(iA, k, h)
                    bsm, auxm = displaced_bsets(iA, k, -h)
                    Pmn_p = ERI_2e3c(bsp, auxp); Jinv_p = inv(ERI_2e2c(auxp))
                    Pmn_m = ERI_2e3c(bsm, auxm); Jinv_m = inv(ERI_2e2c(auxm))
                    @tensor ERI_p[m, n, r, s] := Pmn_p[m, n, A] * Jinv_p[A, B] * Pmn_p[r, s, B]
                    @tensor ERI_m[m, n, r, s] := Pmn_m[m, n, A] * Jinv_m[A, B] * Pmn_m[r, s, B]
                    dERI = (ERI_p .- ERI_m) ./ (2h / B2A)

                    @tensor Jq_ref[m, n] := D[r, s] * dERI[m, n, r, s]
                    @tensor Kq_ref[m, n] := D[r, s] * dERI[m, r, n, s]

                    @test Jq[:, :, k] ≈ Jq_ref atol=1e-5
                    @test Kq[:, :, k] ≈ Kq_ref atol=1e-5
                end
            end
        end

        @testset "CPHF orbital response (finite difference of the SCF density)" begin
            C = wf.orbitals.C
            ndocc = wf.ndocc
            Co = C[:, 1:ndocc]
            Cv = C[:, ndocc+1:end]
            bset = cache.bset

            function displaced_molstring(iA, k, delta)
                newatoms = deepcopy(atoms)
                xyz = collect(newatoms[iA].xyz)
                xyz[k] += delta
                newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
                return join(["$(Int(a.Z))   $(a.xyz[1])   $(a.xyz[2])   $(a.xyz[3])" for a in newatoms], "\n")
            end

            function scf_D(molstring)
                Fermi.Options.set("molstring", molstring)
                Fermi.Options.set("unit", "angstrom")
                Fermi.Options.set("basis", basis)
                Fermi.Options.set("df", true)
                ints_ = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.RIFIT())
                wf_ = Fermi.HartreeFock.RHF(ints_, Fermi.HartreeFock.get_rhf_alg())
                Co_ = wf_.orbitals.C[:, 1:wf_.ndocc]
                return Co_ * Co_'
            end

            B2A = Molecules.bohr_to_angstrom
            h = 5e-4

            for iA in 1:natm
                U = Fermi.HartreeFock.cphf_solve_df(wf, aoints, iA, cache)

                ∂S = zeros(size(Co, 1), size(Co, 1), 3)
                GaussianBasis.∇overlap!(∂S, bset, iA)

                for q in 1:3
                    Dp = scf_D(displaced_molstring(iA, q, h))
                    Dm = scf_D(displaced_molstring(iA, q, -h))
                    dD_FD = (Dp .- Dm) ./ (2h / B2A)

                    ΔD = Cv * U[:, :, q] * Co'
                    ΔD .+= ΔD'
                    S_oo = Co' * (@view ∂S[:, :, q]) * Co
                    dD_analytic = ΔD .- Co * S_oo * Co'

                    @test dD_analytic ≈ dD_FD atol=1e-5
                end
            end

            Fermi.Options.set("molstring", """
            O   0.000000000000   0.000000000000  -0.143225816552
            H   0.000000000000   1.638036840407   1.136548822547
            H   0.000000000000  -1.638036840407   1.136548822547
            """)
            Fermi.Options.set("unit", "bohr")
            Fermi.Options.set("basis", basis)
            Fermi.Options.set("df", true)
        end

        @reset
        @set printstyle none
    end

    @testset "DF direct two-electron Hessian (ERI_hess_JK_df)" begin
        # Phase 4 of the DF-RHF analytic Hessian plan: the DF analog of
        # ERI_hess_JK, i.e. d2X^A_DF/dAdB with P/C held fixed. See DFHess.jl's
        # header comment for the derivation (a from-scratch product-rule
        # rederivation of d2c/dAdB, d2d/dAdB, d2W/dAdB, d2V/dAdB -- an
        # earlier hand-simplified attempt using shortcuts analogous to how
        # the *first* derivative collapses was wrong in two independent
        # ways, caught by this exact checkpoint).
        Fermi.Options.set("molstring", """
        O   0.000000000000   0.000000000000  -0.143225816552
        H   0.000000000000   1.638036840407   1.136548822547
        H   0.000000000000  -1.638036840407   1.136548822547
        """)
        Fermi.Options.set("unit", "bohr")
        basis = "sto-3g"
        Fermi.Options.set("basis", basis)
        Fermi.Options.set("df", true)

        aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.RIFIT())
        wf = Fermi.HartreeFock.RHF(aoints, Fermi.HartreeFock.get_rhf_alg())
        cache = Fermi.HartreeFock.build_df_hess_cache(wf, aoints)
        atoms = wf.molecule.atoms
        natm = length(atoms)

        ndocc = wf.ndocc
        Co = wf.orbitals.C[:, 1:ndocc]
        P = 2.0 .* (Co * Co')

        @testset "Block symmetry (free, FD-independent)" begin
            for iA in 1:natm, iB in 1:natm
                HAB = Fermi.HartreeFock.ERI_hess_JK_df(cache, P, iA, iB)
                HBA = Fermi.HartreeFock.ERI_hess_JK_df(cache, P, iB, iA)
                @test HAB ≈ HBA' atol=1e-8
            end
        end

        @testset "vs double finite difference of the dense DF direct energy" begin
            function X_DF(newatoms)
                bs = BasisSet(basis, newatoms)
                aux = BasisSet(aoints.eri_type.basisset.name, newatoms)
                Pmn_ = ERI_2e3c(bs, aux)
                Jinv_ = inv(ERI_2e2c(aux))
                @tensor ERI[m, n, r, s] := Pmn_[m, n, A] * Jinv_[A, B] * Pmn_[r, s, B]
                @tensor Xj = 0.5 * P[m, n] * P[r, s] * ERI[m, n, r, s]
                @tensor Xk = 0.25 * P[m, n] * P[r, s] * ERI[m, s, r, n]
                return Xj - Xk
            end

            function displaced(atoms, iA, kA, dA, iB, kB, dB)
                newatoms = deepcopy(atoms)
                xyzA = collect(newatoms[iA].xyz); xyzA[kA] += dA
                newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyzA)
                xyzB = collect(newatoms[iB].xyz); xyzB[kB] += dB
                newatoms[iB] = Molecules.Atom(newatoms[iB].Z, newatoms[iB].mass, xyzB)
                return newatoms
            end

            h = 5e-3
            B2A = Molecules.bohr_to_angstrom

            for iA in 1:natm, iB in 1:natm
                Hd = Fermi.HartreeFock.ERI_hess_JK_df(cache, P, iA, iB)
                for qA in 1:3, qB in 1:3
                    if iA == iB && qA == qB
                        fp = X_DF(displaced(atoms, iA, qA, h, iB, qB, 0.0))
                        f0 = X_DF(atoms)
                        fm = X_DF(displaced(atoms, iA, qA, -h, iB, qB, 0.0))
                        d2 = (fp - 2 * f0 + fm) / (h / B2A)^2
                    else
                        fpp = X_DF(displaced(atoms, iA, qA, h, iB, qB, h))
                        fpm = X_DF(displaced(atoms, iA, qA, h, iB, qB, -h))
                        fmp = X_DF(displaced(atoms, iA, qA, -h, iB, qB, h))
                        fmm = X_DF(displaced(atoms, iA, qA, -h, iB, qB, -h))
                        d2 = (fpp - fpm - fmp + fmm) / (4 * (h / B2A)^2)
                    end
                    @test Hd[qA, qB] ≈ d2 atol=1e-3
                end
            end
        end

        @reset
        @set printstyle none
    end

    @testset "Full DF-RHF Hessian (finite difference of the DF gradient, and Psi4)" begin
        # Phase 5 of the DF-RHF analytic Hessian plan: end-to-end assembly
        # (RHFhess dispatching on an AbstractDFERI aoints), combining the
        # DF direct term (ERI_hess_JK_df) with CPHF orbital/density response
        # built from cphf_solve_df/eri_grad_JK_df.
        #
        # Aux basis forced to cc-pvtz-jkfit rather than RIFIT's "auto"
        # fallback (cc-pvqz-rifit for a non-Dunning primary basis like
        # sto-3g): that fallback is a poor size match for a 7-function
        # primary basis (242 aux functions), giving a metric condition
        # number ~7.7e6 that's fine for a first derivative (DFgrad.jl) but
        # amplifies ordinary cross-run floating-point/threading noise (CG
        # iteration order, etc.) past 1e-3 once squared through a *second*
        # derivative and an independently-rebuilt aoints/CPHF solve --
        # caught by this test itself (both sub-checks below originally
        # failed at the ~1e-3 level with the auto aux basis, vanishing once
        # switched to this better-conditioned one).
        path = joinpath(@__DIR__, "xyz/water.xyz")
        mol = open(f->read(f,String), path)
        basis = "sto-3g"

        Fermi.Options.set("molstring", mol)
        Fermi.Options.set("basis", basis)
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvtz-jkfit")

        aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.JKFIT())
        wf = Fermi.HartreeFock.RHF(aoints, Fermi.HartreeFock.get_rhf_alg())
        H = Fermi.HartreeFock.RHFhess(aoints, wf)

        @test H ≈ H' atol=1e-8

        # convenience dispatch (RHFhess(wf), no aoints) must pick the same
        # DF path automatically via the `df` option -- IntegralHelper's own
        # default eri_type selection (JKFIT/RIFIT vs Chonky/SparseERI).
        H_conv = Fermi.HartreeFock.RHFhess(wf)
        @test H ≈ H_conv atol=1e-8

        # Independent check: finite difference of Fermi's own DF analytic
        # gradient (DFgrad.jl). Mirrors the exact-ERI "Full RHF Hessian"
        # test's own FD-of-gradient checkpoint exactly (same water.xyz,
        # same default angstrom convention -- atoms are re-stringified from
        # `atoms[iA].xyz`, which Molecules.jl always stores internally in
        # Angstrom regardless of the original input unit, so no unit needs
        # to be re-set here).
        function displaced_molstring(atoms, iA, k, delta)
            newatoms = deepcopy(atoms)
            xyz = collect(newatoms[iA].xyz)
            xyz[k] += delta
            newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
            return join(["$(Int(a.Z))   $(a.xyz[1])   $(a.xyz[2])   $(a.xyz[3])" for a in newatoms], "\n")
        end

        function scf_grad_df(molstring)
            Fermi.Options.set("molstring", molstring)
            Fermi.Options.set("basis", basis)
            Fermi.Options.set("df", true)
            Fermi.Options.set("jkfit", "cc-pvtz-jkfit")
            ints_ = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.JKFIT())
            wf_ = Fermi.HartreeFock.RHF(ints_, Fermi.HartreeFock.get_rhf_alg())
            return Fermi.HartreeFock.RHFgrad(ints_, wf_)
        end

        atoms = wf.molecule.atoms
        natm = length(atoms)
        B2A = Molecules.bohr_to_angstrom
        h = 5e-4
        H_FD = zeros(3*natm, 3*natm)
        for iB in 1:natm, qB in 1:3
            gp = scf_grad_df(displaced_molstring(atoms, iB, qB, h))
            gm = scf_grad_df(displaced_molstring(atoms, iB, qB, -h))
            H_FD[:, 3*(iB-1)+qB] .= vec(((gp .- gm) ./ (2h/B2A))')
        end
        @test H ≈ H_FD atol=1e-4

        @reset
        @set printstyle none
    end

    @testset "Full DF-RHF Hessian vs Psi4" begin
        # Independent-implementation cross-check (not just finite
        # difference): Psi4 1.10 `hessian('scf')`, scf_type df, water/sto-3g,
        # noreorient/nocom so atom/axis ordering matches, this exact
        # geometry (bohr). Aux basis must match on both sides to be
        # meaningful -- neither side's "auto" heuristic picks the same
        # basis (Fermi's RIFIT "auto" defaults to cc-pvqz-rifit for a
        # non-Dunning primary basis; Psi4's own DF_BASIS_SCF default is a
        # different heuristic entirely), and Psi4's analytic Hessian engine
        # rejects aux bases with angular momentum > g (cc-pvqz-jkfit's h
        # functions raise `lmax_exceeded`) -- so cc-pvtz-jkfit (max L=4) is
        # forced explicitly on both sides via the `jkfit` option / Psi4's
        # `df_basis_scf`.
        Fermi.Options.set("molstring", """
        O   0.000000000000   0.000000000000  -0.143225816552
        H   0.000000000000   1.638036840407   1.136548822547
        H   0.000000000000  -1.638036840407   1.136548822547
        """)
        Fermi.Options.set("unit", "bohr")
        Fermi.Options.set("basis", "sto-3g")
        Fermi.Options.set("df", true)
        Fermi.Options.set("jkfit", "cc-pvtz-jkfit")
        Fermi.Options.set("scf_e_conv", 1e-12)
        Fermi.Options.set("scf_max_rms", 1e-10)
        Fermi.Options.set("cphf_conv", 1e-10)

        aoints = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.JKFIT())
        wf = Fermi.HartreeFock.RHF(aoints, Fermi.HartreeFock.get_rhf_alg())
        H = Fermi.HartreeFock.RHFhess(aoints, wf)

        @test H ≈ H' atol=1e-9

        H_psi4 = [
            0.07613774603504    -0.00000000000000    -0.00000000000000    -0.03806887301680     0.00000000000000     0.00000000000000    -0.03806887301707    -0.00000000000000    -0.00000000000000
           -0.00000000000000     0.48291234123351    -0.00000000000038     0.00000000000000    -0.24145617061195    -0.15890357222942    -0.00000000000000    -0.24145617062065     0.15890357222497
           -0.00000000000000    -0.00000000000038     0.43735805543467     0.00000000000000    -0.07344719886939    -0.21867902771937    -0.00000000000000     0.07344719887495    -0.21867902771719
           -0.03806887301680     0.00000000000000     0.00000000000000     0.04537590838802    -0.00000000000000    -0.00000000000000    -0.00730703537143    -0.00000000000000     0.00000000000000
            0.00000000000000    -0.24145617061195    -0.07344719886939    -0.00000000000000     0.25786921957443     0.11617538555496    -0.00000000000000    -0.01641304895170    -0.04272818667789
            0.00000000000000    -0.15890357222942    -0.21867902771937    -0.00000000000000     0.11617538555496     0.19775552422997    -0.00000000000000     0.04272818667623     0.02092350348770
           -0.03806887301707    -0.00000000000000    -0.00000000000000    -0.00730703537143    -0.00000000000000    -0.00000000000000     0.04537590838906     0.00000000000000     0.00000000000000
           -0.00000000000000    -0.24145617062065     0.07344719887495    -0.00000000000000    -0.01641304895170     0.04272818667623     0.00000000000000     0.25786921956403    -0.11617538554872
           -0.00000000000000     0.15890357222497    -0.21867902771719     0.00000000000000    -0.04272818667789     0.02092350348770     0.00000000000000    -0.11617538554872     0.19775552423506
        ]
        @test H ≈ H_psi4 atol=1e-5

        @reset
        @set printstyle none
    end

    @testset "Full RHF Hessian (finite difference of the gradient, and Psi4)" begin
        Fermi.Options.set("molstring", """
        O   0.000000000000   0.000000000000   0.000000000000
        H   0.758602190000   0.000000000000   0.504284980000
        H   0.758602190000   0.000000000000  -0.504284980000
        """)
        basis = "sto-3g"
        Fermi.Options.set("basis", basis)
        wf = @energy rhf
        H = Fermi.HartreeFock.RHFhess(wf)

        @test H ≈ H' atol=1e-10

        # Reference from Psi4 1.10 (`hessian('scf')`, basis sto-3g, this exact
        # geometry, noreorient/nocom so the atom/axis ordering matches).
        H_psi4 = [
            1.13060065634761   0.0                0.0               -0.56530032817424   0.0               -0.44471126604731  -0.56530032817424   0.0                0.44471126604731
            0.0               -0.07707295928527   0.0                0.0                0.03853647964198   0.0                0.0                0.03853647964198   0.0
            0.0                0.0                0.61840092098178  -0.52310501956430   0.0               -0.30920046049119   0.52310501956430   0.0               -0.30920046049119
           -0.56530032817424   0.0               -0.52310501956430   0.56393084645098   0.0                0.48390814280581   0.00136948172326   0.0                0.03919687675849
            0.0                0.03853647964198   0.0                0.0               -0.13144426267318   0.0                0.0                0.09290778303120   0.0
           -0.44471126604731   0.0               -0.30920046049119   0.48390814280581   0.0                0.53179651540804  -0.03919687675849   0.0               -0.22259605491685
           -0.56530032817424   0.0                0.52310501956430   0.00136948172326   0.0               -0.03919687675849   0.56393084645098   0.0               -0.48390814280581
            0.0                0.03853647964198   0.0                0.0                0.09290778303120   0.0                0.0               -0.13144426267318   0.0
            0.44471126604731   0.0               -0.30920046049119   0.03919687675849   0.0               -0.22259605491685  -0.48390814280581   0.0                0.53179651540804
        ]
        @test H ≈ H_psi4 atol=1e-6

        # Independent check: finite difference of Fermi's own analytic gradient.
        function displaced_molstring(atoms, iA, k, delta)
            newatoms = deepcopy(atoms)
            xyz = collect(newatoms[iA].xyz)
            xyz[k] += delta
            newatoms[iA] = Molecules.Atom(newatoms[iA].Z, newatoms[iA].mass, xyz)
            return join(["$(Int(a.Z))   $(a.xyz[1])   $(a.xyz[2])   $(a.xyz[3])" for a in newatoms], "\n")
        end

        function scf_grad(molstring)
            Fermi.Options.set("molstring", molstring)
            Fermi.Options.set("basis", basis)
            wf_ = @energy rhf
            return Fermi.HartreeFock.RHFgrad(wf_)
        end

        atoms = wf.molecule.atoms
        natm = length(atoms)
        B2A = Molecules.bohr_to_angstrom
        h = 5e-4
        H_FD = zeros(3*natm, 3*natm)
        for iB in 1:natm, qB in 1:3
            gp = scf_grad(displaced_molstring(atoms, iB, qB, h))
            gm = scf_grad(displaced_molstring(atoms, iB, qB, -h))
            H_FD[:, 3*(iB-1)+qB] .= vec(((gp .- gm) ./ (2h/B2A))')
        end
        @test H ≈ H_FD atol=1e-5
    end

    @testset "@hessian macro and vibrational_analysis" begin
        # Harmonic vibrational analysis via simple mass-weighted projection is
        # only rigorously defined AT a stationary point (zero gradient) --
        # away from one, translation/rotation contaminate the Hessian in a
        # way projection doesn't fully remove, and even Psi4's own frequency
        # routine visibly breaks (it explicitly warns "Vibrations include
        # un-projected rotation-like modes" and reports spurious large
        # imaginary frequencies) if you feed it a non-equilibrium geometry.
        # So: optimize first, matching what any real use of this would do.
        Fermi.Options.set("molstring", """
        O   0.000000000000   0.000000000000   0.000000000000
        H   0.758602190000   0.000000000000   0.504284980000
        H   0.758602190000   0.000000000000  -0.504284980000
        """)
        Fermi.Options.set("basis", "sto-3g")
        Fermi.Options.set("geom_e_conv", 1e-10)
        Fermi.Options.set("geom_grms_conv", 1e-9)
        Fermi.Options.set("geom_gmax_conv", 1e-9)
        Fermi.Options.set("geom_drms_conv", 1e-8)
        Fermi.Options.set("geom_dmax_conv", 1e-8)
        Fermi.Options.set("geom_max_iter", 100)
        wf = @optimize rhf
        @test maximum(abs.(Fermi.HartreeFock.RHFgrad(wf))) < 1e-8

        H = @hessian wf => rhf
        freqs, modes = vibrational_analysis(H, wf.molecule)

        # 6 translation/rotation modes (nonlinear triatomic) should be
        # essentially zero; the remaining 3 are the genuine vibrations.
        @test count(f -> abs(f) < 1.0, freqs) == 6

        vib_freqs = sort(filter(f -> abs(f) >= 1.0, freqs))
        # Reference: Psi4 1.10 `frequencies('scf')`, sto-3g, at this exact
        # (Fermi-optimized) geometry -- Psi4's own TR/V labeling comes out
        # clean here (unlike at a non-stationary geometry), confirming this
        # is the right kind of test. atol=1.0 cm^-1 comfortably covers the
        # ~0.2-0.4 cm^-1 residual expected from ordinary cross-code numerical
        # noise (integral thresholds, summation order, etc.) at this level of
        # agreement -- already essentially machine-converged (gradient <1e-8).
        @test vib_freqs ≈ [2170.04597272, 4140.00188775, 4391.06663761] atol=1.0

        # @frequencies combines the two calls above into one; @frequency is
        # a plain alias for it (same behavior, different name).
        freqs_f, modes_f = @frequencies wf => rhf
        @test freqs_f ≈ freqs atol=1e-10
        @test modes_f ≈ modes atol=1e-10

        freqs_f2, modes_f2 = @frequency wf => rhf
        @test freqs_f2 ≈ freqs atol=1e-10
        @test modes_f2 ≈ modes atol=1e-10

        # Linear molecule sanity check: only 5 (not 6) trans/rot zero modes,
        # since rotation about the molecular axis isn't a physical DOF.
        Fermi.Options.set("molstring", """
        N   0.000000000000   0.000000000000   0.000000000000
        N   0.000000000000   0.000000000000   1.100000000000
        """)
        wf_n2 = @energy rhf
        H_n2 = @hessian wf_n2 => rhf
        freqs_n2, _ = vibrational_analysis(H_n2, wf_n2.molecule)
        @test count(f -> abs(f) < 1.0, freqs_n2) == 5
    end
end

@reset
