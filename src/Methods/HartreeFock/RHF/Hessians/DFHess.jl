using LinearAlgebra
using TensorOperations

# DF-RHF analytic Hessian. Extends CPHF.jl/RHFhess.jl the same way DFgrad.jl
# extends RHFgrad.jl -- see DFgrad.jl's header for the underlying RI-HF
# energy expression and the already-validated (finite difference, ~1e-8,
# first attempt) Coulomb/exchange gradient formulas this builds on:
#
#   c_P := D_mn*(mn|P), d_P := J^-1_PQ*c_Q            (Coulomb)
#   W[m,i,P] := sum_n C_ni*(mn|P), V := W*J^-1         (exchange, i=occupied)
#
# Phase 3 (this section): Jq_DF/Kq_DF, the DF analog of CPHF.jl's
# eri_grad_JK -- matrix-valued (not summed against a second D/P) first
# derivatives, needed for cphf_rhs's "skeleton" term. Derivation (worked out
# and independently re-derived during planning, see the DF-RHF Analytic
# Hessian plan's Theory section):
#
#   Jq_DF[m,n] = sum_P [d(mn|P)/dA]*d_P + sum_P (mn|P)*(dd_P/dA - f_P)
#   Kq_DF[m,n] = sum_i,P dW[m,i,P]*V[n,i,P] + sum_i,P W[m,i,P]*dV_expl[n,i,P]
#                - sum_i V_i[m,:]'*(dJ/dA)*V_i[n,:]
#
# where dc_P/dA := sum_mn D_mn*d(mn|P)/dA, dd_P/dA := J^-1*(dc/dA),
# f_P := J^-1*(dJ/dA)*d, dW[m,i,P] := sum_n C_ni*d(mn|P)/dA (same C, fixed --
# no CPHF/orbital response needed for a FIRST derivative, same reasoning as
# the gradient), dV_expl := dW*J^-1 (the "explicit" part of dV/dA, not
# including the J^-1-derivative piece, which the last Kq_DF term handles
# separately -- exactly mirroring how DFgrad.jl splits the Coulomb term into
# dc.d and -0.5*d'(dJ/dA)d).
#
# Free, FD-independent correctness check: Kq_DF must come out symmetric
# (Kq_DF ≈ Kq_DF'), same as the exact-ERI Kq is provably symmetric --
# swapping m<->n in the formula above swaps the first two terms into each
# other and leaves the third manifestly symmetric (dJ/dA is a symmetric
# matrix). Checked as part of this phase's validation.
#
# eri_grad_JK (CPHF.jl) itself is untouched -- cphf_rhs/_rhf_hess_response
# (CPHF.jl/RHFhess.jl) instead get DF-specific counterparts here
# (cphf_rhs_df/cphf_solve_df) that call eri_grad_JK_df, reusing
# cphf_Amatvec/build_response_fock! unchanged (they only touch build_fock!,
# already eri_type-generic).

"""
    build_df_hess_cache(wfn::RHF, ints::IntegralHelper{Float64,<:AbstractDFERI})

Atom-independent quantities for the DF-RHF Hessian, built once and reused
across every atom: the 3-center/2-center integral tensors (`Pmn`, `Jmat`),
the metric inverse (`Jinv`), and the occupied-orbital half-transform (`W`,
`V`) -- exactly the objects `DFgrad.jl` builds before its own atom loop,
reused here to avoid recomputing an `naux x naux` inversion or the `W`
half-transform once per atom.
"""
function build_df_hess_cache(wfn::RHF, ints::Fermi.Integrals.IntegralHelper{Float64,<:Fermi.Integrals.AbstractDFERI})
    bset = BasisSet(ints.basis, ints.molecule.atoms)
    auxbset = ints.eri_type.basisset
    ndocc = wfn.ndocc
    Co = wfn.orbitals.C[:, 1:ndocc]
    D = Co * Co'

    Pmn = GaussianBasis.ERI_2e3c(bset, auxbset)
    Jmat = GaussianBasis.ERI_2e2c(auxbset)
    Jinv = inv(Jmat)

    @tensor c[A] := D[m, n] * Pmn[m, n, A]
    d = Jinv * c

    @tensor W[m, i, A] := Co[n, i] * Pmn[m, n, A]
    @tensor V[m, i, A] := W[m, i, B] * Jinv[B, A]

    # Qtot_PQ := sum_i (D*V_i)'*V_i, a small (naux,naux) matrix. Atom-
    # independent (built entirely from the converged, zeroth-order D/V), so
    # it's computed once here and reused as-is by Phase 4's direct term --
    # differentiating Tr(Qtot.dJ/dA) w.r.t. a second atom only ever hits
    # dJ/dA, never Qtot itself.
    naux = auxbset.nbas
    Qtot = zeros(naux, naux)
    for i in 1:ndocc
        Vi = @view V[:, i, :]
        Qtot .+= (D * Vi)' * Vi
    end

    return (bset=bset, auxbset=auxbset, Co=Co, D=D, Pmn=Pmn, Jmat=Jmat, Jinv=Jinv,
            c=c, d=d, W=W, V=V, Qtot=Qtot)
end

"""
    df_hess_atom_derivs(cache, iA)

Per-atom first-derivative quantities for atom `iA`, all three Cartesian
directions: `Jq`/`Kq` (see this file's header), plus the intermediates
(`dW`, `dV_expl`, `∇Jmat`) Phase 4's direct second-derivative term reuses
rather than recomputing. `eri_grad_JK_df` is a thin wrapper around this
returning just `(Jq, Kq)`.
"""
function df_hess_atom_derivs(cache, iA::Int)
    bset, auxbset = cache.bset, cache.auxbset
    D, Co = cache.D, cache.Co
    Pmn, Jinv, d, W, V = cache.Pmn, cache.Jinv, cache.d, cache.W, cache.V
    nbas = bset.nbas
    naux = auxbset.nbas

    ∇Pmn = zeros(nbas, nbas, naux, 3)
    GaussianBasis.∇ERI_2e3c!(∇Pmn, bset, auxbset, iA)
    ∇Jmat = zeros(naux, naux, 3)
    GaussianBasis.∇ERI_2e2c!(∇Jmat, auxbset, iA)

    Jq = zeros(nbas, nbas, 3)
    Kq = zeros(nbas, nbas, 3)
    dW_all = zeros(nbas, size(Co, 2), naux, 3)
    dV_expl_all = zeros(nbas, size(Co, 2), naux, 3)
    dc_all = zeros(naux, 3)
    dd_all = zeros(naux, 3)

    for q in 1:3
        ∇Pq = @view ∇Pmn[:, :, :, q]
        ∇Jq = @view ∇Jmat[:, :, q]

        @tensor dc[A] := D[m, n] * ∇Pq[m, n, A]
        dd = Jinv * dc
        f = Jinv * (∇Jq * d)
        dc_all[:, q] .= dc
        dd_all[:, q] .= dd

        @tensor Jq_q[m, n] := ∇Pq[m, n, A] * d[A] + Pmn[m, n, A] * (dd[A] - f[A])
        Jq[:, :, q] .= Jq_q

        @tensor dW[m, i, A] := Co[n, i] * ∇Pq[m, n, A]
        @tensor dV_expl[m, i, A] := dW[m, i, B] * Jinv[B, A]
        dW_all[:, :, :, q] .= dW
        dV_expl_all[:, :, :, q] .= dV_expl

        @tensor Kq_q[m, n] := dW[m, i, A] * V[n, i, A] + W[m, i, A] * dV_expl[n, i, A]
        @tensor Kcorr[m, n] := V[m, i, A] * ∇Jq[A, B] * V[n, i, B]
        Kq[:, :, q] .= Kq_q .- Kcorr
    end

    return (Jq=Jq, Kq=Kq, ∇Pmn=∇Pmn, ∇Jmat=∇Jmat, dW=dW_all, dV_expl=dV_expl_all,
            dc=dc_all, dd=dd_all)
end

eri_grad_JK_df(cache, iA::Int) = (r = df_hess_atom_derivs(cache, iA); (r.Jq, r.Kq))

"""
    cphf_rhs_df(wfn, ints, iA, cache)

DF analog of `cphf_rhs` (CPHF.jl) -- identical structure, using
`eri_grad_JK_df`/`cache` instead of `eri_grad_JK`/`ij_vals`/`σvals`.
`build_fock!` (via `ints`) is unchanged and already dispatches on DF.
"""
function cphf_rhs_df(wfn::RHF, ints, iA::Int, cache)
    bset = cache.bset
    nbas = bset.nbas
    ndocc = wfn.ndocc
    nvir = nbas - ndocc
    C = wfn.orbitals.C
    Co = @view C[:, 1:ndocc]
    Cv = @view C[:, ndocc+1:end]
    eps_o = wfn.orbitals.eps[1:ndocc]

    ∂H = zeros(nbas, nbas, 3)
    ∂Vtmp = zeros(nbas, nbas, 3)
    GaussianBasis.∇kinetic!(∂H, bset, iA)
    GaussianBasis.∇nuclear!(∂Vtmp, bset, iA)
    ∂H .+= ∂Vtmp

    ∂S = zeros(nbas, nbas, 3)
    GaussianBasis.∇overlap!(∂S, bset, iA)

    Jq_all, Kq_all = eri_grad_JK_df(cache, iA)

    B = zeros(nvir, ndocc, 3)
    G = zeros(nbas, nbas)
    H0 = zeros(nbas, nbas)
    for q in 1:3
        Fskel = @view(∂H[:, :, q]) .+ 2 .* (@view Jq_all[:, :, q]) .- (@view Kq_all[:, :, q])

        Sq = @view ∂S[:, :, q]
        S_oo = Co' * Sq * Co
        dD_S = -Co * S_oo * Co'
        build_fock!(G, H0, dD_S, ints)

        B[:, :, q] .= Cv' * Fskel * Co .- (Cv' * Sq * Co) * Diagonal(eps_o) .+ Cv' * G * Co
    end

    return B, ∂H, ∂S, Jq_all, Kq_all
end

"""
    cphf_solve_df(wfn, ints, iA, cache)

DF analog of `cphf_solve` (CPHF.jl). Thin wrapper around
`cphf_solve_df_full` for callers that only need `U`.
"""
function cphf_solve_df(wfn::RHF, ints, iA::Int, cache)
    U, _, _, _, _ = cphf_solve_df_full(wfn, ints, iA, cache)
    return U
end

"""
    cphf_solve_df_full(wfn, ints, iA, cache)

DF analog of `cphf_solve_full` (CPHF.jl): solves the Jacobi-preconditioned
CPHF equations and also returns `cphf_rhs_df`'s intermediates alongside
`U`. See `cphf_solve_full`'s comment for the preconditioning derivation
(`cphf_Amatvec` is shared, unchanged, between the exact-ERI and DF paths).
"""
function cphf_solve_df_full(wfn::RHF, ints, iA::Int, cache)
    ndocc = wfn.ndocc
    nbas = size(wfn.orbitals.C, 1)
    nvir = nbas - ndocc
    C = wfn.orbitals.C
    Co = @view C[:, 1:ndocc]
    Cv = @view C[:, ndocc+1:end]
    eps = wfn.orbitals.eps
    Δeps = eps[ndocc+1:end] .- eps[1:ndocc]'
    d = sqrt.(Δeps)

    maxiter = Options.get("cphf_max_iter")
    tol = Options.get("cphf_conv")

    B, ∂H, ∂S, Jq_all, Kq_all = cphf_rhs_df(wfn, ints, iA, cache)
    U = zeros(nvir, ndocc, 3)
    for q in 1:3
        precond_matvec(y) = y .+ cphf_Amatvec(y ./ d, Co, Cv, ints) ./ d
        rhs = -B[:, :, q] ./ d
        sol, info = KrylovKit.linsolve(precond_matvec, rhs, zeros(nvir, ndocc), KrylovKit.CG(; maxiter=maxiter, tol=tol))
        U[:, :, q] .= sol ./ d
    end
    return U, ∂H, ∂S, Jq_all, Kq_all
end

# --- Phase 4: DF "direct" (P frozen) second-derivative term ---
#
# DF analog of ERI_hess_JK (TwoElectronHess.jl): d2X^A_DF/dAdB with P held
# fixed at its converged value. Coulomb uses P-based c,d (matching
# DFgrad.jl's ORIGINAL gradient convention -- 0.5*c^P.d^P is *directly*
# X^A's 0.5*P.P.(mn|rs) contribution); since P=2D, c^P=2c^D, d^P=2d^D
# exactly, reusing df_hess_atom_derivs's D-based dc via a factor of 2.
#
# An earlier version of this derivation tried to hand-simplify the second
# derivative into a minimal closed form (mirroring how the *first*
# derivative collapses via m<->n relabeling into the compact Jq_DF/Kq_DF
# expressions) and got it wrong twice: it treated Phase 3's cached "dd"
# (Jinv*dc only, deliberately missing the metric-derivative piece so
# Jq_DF's own formula could add it back separately) as if it were the FULL
# derivative of d, and treated Qtot (built from the converged, zeroth-order
# V) as if it had no further B-dependence -- but V=W*Jinv itself depends on
# geometry, so Qtot=sum_i V_i'DV_i does too, and differentiating
# Tr(Qtot.dJ/dA) w.r.t. B needs a dQtot/dB term that formula silently
# dropped. Caught by this phase's own finite-difference checkpoint (a ~4
# Hartree/bohr^2 discrepancy on water/sto-3g, obviously wrong -- not a
# subtle bug). Fixed by computing every second derivative object directly
# via the product/chain rule instead of hand-simplifying, treating d and V
# each as ordinary functions of geometry with their own (already-validated
# at first order) derivative rules:
#
#   dJinv(A) = -Jinv.∇J_A.Jinv
#   d2Jinv(A,B) = Jinv.∇J_A.Jinv.∇J_B.Jinv + Jinv.∇J_B.Jinv.∇J_A.Jinv - Jinv.∇2J_AB.Jinv
#
#   dd(A)_full = Jinv.dc(A) + dJinv(A).c            (Phase 3's cached "dd" is only the first piece)
#   d2d(A,B) = Jinv.d2c(A,B) + dJinv(A).dc(B) + dJinv(B).dc(A) + d2Jinv(A,B).c
#   E_J = 0.5*c.d  =>  d2E_J/dAdB = 0.5*[d2c(A,B).d + dc(A).dd(B)_full + dc(B).dd(A)_full + c.d2d(A,B)]
#
#   dV(A) = dW(A).Jinv + W.dJinv(A) = dV_expl(A) - g(A), g(A) := V.∇J_A.Jinv  (Phase 3)
#   d2V(A,B) = d2W(A,B).Jinv + dW(A).dJinv(B) + dW(B).dJinv(A) + W.d2Jinv(A,B)
#   E_K = -sum_i D.W.V  =>  d2E_K/dAdB = -sum_i D.[d2W(A,B).V + dW(A).dV(B) + dW(B).dV(A) + W.d2V(A,B)]
#
# (both E_J's and E_K's *first*-derivative forms were re-derived this same
# way as a self-check and confirmed to reduce to the already-validated
# gradient formulas before trusting the second-derivative extension).
#
# d2W[m,i,P](A,B) := sum_n C[n,i]*d2(mn|P)/dAdB -- not a new kernel, just
# the new ∇2ERI_2e3c Hessian integral half-transformed by the converged
# (fixed) C, computed fresh per atom pair (not cached -- an
# O(natm^2)-sized cache of an (nbas,ndocc,naux) object per pair would be
# wasteful, mirroring how ERI_hess_JK itself works fresh per pair).
#
# X^A_DF's second derivative is d2E_J/dAdB + d2E_K/dAdB directly (same
# +dE_K/dA relationship established for the gradient/Phase 3, no extra
# sign or prefactor -- K's own -4/-0.25 factors already cancel exactly).

"""
    ERI_hess_JK_df(cache, P::AbstractMatrix, iA::Int, iB::Int)

DF analog of `ERI_hess_JK` (`TwoElectronHess.jl`) -- the two-electron
Coulomb+exchange "direct" (density held fixed) contribution to the DF-RHF
Hessian block `(iA,iB)`. See this file's header comment above for the full
derivation. `P` is the converged density (`2*Co*Co'`); `cache` is
`build_df_hess_cache`'s output.
"""
function ERI_hess_JK_df(cache, P::AbstractMatrix, iA::Int, iB::Int)
    bset, auxbset = cache.bset, cache.auxbset
    D, Co, Jinv, W, V = cache.D, cache.Co, cache.Jinv, cache.W, cache.V
    cP = 2.0 .* cache.c
    dP = 2.0 .* cache.d

    rA = df_hess_atom_derivs(cache, iA)
    rB = df_hess_atom_derivs(cache, iB)

    ∇2Pmn = zeros(bset.nbas, bset.nbas, auxbset.nbas, 3, 3)
    GaussianBasis.∇2ERI_2e3c!(∇2Pmn, bset, auxbset, iA, iB)
    ∇2Jmat = zeros(auxbset.nbas, auxbset.nbas, 3, 3)
    GaussianBasis.∇2ERI_2e2c!(∇2Jmat, auxbset, iA, iB)

    H = zeros(3, 3)
    for qA in 1:3, qB in 1:3
        ∇Jq_A = @view rA.∇Jmat[:, :, qA]
        ∇Jq_B = @view rB.∇Jmat[:, :, qB]
        ∇2Pq = @view ∇2Pmn[:, :, :, qA, qB]
        ∇2Jq = @view ∇2Jmat[:, :, qA, qB]

        dJinv_A = -Jinv * ∇Jq_A * Jinv
        dJinv_B = -Jinv * ∇Jq_B * Jinv
        d2Jinv_AB = Jinv * ∇Jq_A * Jinv * ∇Jq_B * Jinv .+ Jinv * ∇Jq_B * Jinv * ∇Jq_A * Jinv .-
                    Jinv * ∇2Jq * Jinv

        ### Coulomb
        dcP_A = 2.0 .* (@view rA.dc[:, qA])
        dcP_B = 2.0 .* (@view rB.dc[:, qB])
        ddP_A = Jinv * dcP_A .+ dJinv_A * cP
        ddP_B = Jinv * dcP_B .+ dJinv_B * cP

        @tensor d2cP[A] := P[m, n] * ∇2Pq[m, n, A]
        d2dP = Jinv * d2cP .+ dJinv_A * dcP_B .+ dJinv_B * dcP_A .+ d2Jinv_AB * cP

        Jpart = 0.5 * (dot(d2cP, dP) + dot(dcP_A, ddP_B) + dot(dcP_B, ddP_A) + dot(cP, d2dP))

        ### Exchange
        dW_A = @view rA.dW[:, :, :, qA]
        dW_B = @view rB.dW[:, :, :, qB]

        @tensor g_A[m, i, P] := V[m, i, R] * ∇Jq_A[R, S] * Jinv[S, P]
        @tensor g_B[m, i, P] := V[m, i, R] * ∇Jq_B[R, S] * Jinv[S, P]
        dV_A = (@view rA.dV_expl[:, :, :, qA]) .- g_A
        dV_B = (@view rB.dV_expl[:, :, :, qB]) .- g_B

        @tensor d2W[m, i, P] := Co[n, i] * ∇2Pq[m, n, P]
        @tensor d2V_1[m, i, P] := d2W[m, i, Q] * Jinv[Q, P]
        @tensor d2V_2[m, i, P] := dW_A[m, i, Q] * dJinv_B[Q, P]
        @tensor d2V_3[m, i, P] := dW_B[m, i, Q] * dJinv_A[Q, P]
        @tensor d2V_4[m, i, P] := W[m, i, Q] * d2Jinv_AB[Q, P]
        d2V = d2V_1 .+ d2V_2 .+ d2V_3 .+ d2V_4

        @tensor term1 = D[m, n] * d2W[m, i, P] * V[n, i, P]
        @tensor term2 = D[m, n] * dW_A[m, i, P] * dV_B[n, i, P]
        @tensor term3 = D[m, n] * dW_B[m, i, P] * dV_A[n, i, P]
        @tensor term4 = D[m, n] * W[m, i, P] * d2V[n, i, P]

        Kpart = -(term1 + term2 + term3 + term4)

        H[qA, qB] = Jpart + Kpart
    end

    return H
end

# --- Phase 5: assembly / top-level dispatch ---
#
# DF analog of RHFhess.jl's main assembly loop -- same formula, same
# RHFHessResponse struct (defined in RHFhess.jl, included before this file),
# just built from cphf_solve_df/eri_grad_JK_df instead of
# cphf_solve/eri_grad_JK, and ERI_hess_JK_df instead of ERI_hess_JK for the
# direct two-electron term. build_fock! (CPHF's own Fock builds, and the
# skeleton-density G build below) already dispatches on `ints`'s eri_type,
# so no changes needed there.

function _rhf_hess_response_df(wfn::RHF, ints, cache, iA, Co, Cv)
    bset = cache.bset
    nbas = bset.nbas
    U, ∂H, ∂S, Jq_all, Kq_all = cphf_solve_df_full(wfn, ints, iA, cache)

    dD = zeros(nbas, nbas, 3)
    dF = zeros(nbas, nbas, 3)
    H0 = zeros(nbas, nbas)
    for q in 1:3
        Fskel = @view(∂H[:, :, q]) .+ 2 .* (@view Jq_all[:, :, q]) .- (@view Kq_all[:, :, q])

        Sq = @view ∂S[:, :, q]
        S_oo = Co' * Sq * Co
        dD_S = -Co * S_oo * Co'
        G = zeros(nbas, nbas)
        build_fock!(G, H0, dD_S, ints)

        ΔF = zeros(nbas, nbas)
        build_response_fock!(ΔF, U[:, :, q], Co, Cv, ints)

        dDq = Cv * U[:, :, q] * Co'
        dDq .+= dDq'
        dDq .+= dD_S
        dD[:, :, q] .= dDq
        dF[:, :, q] .= Fskel .+ ΔF .+ G
    end

    return RHFHessResponse(dD, dF, ∂H, ∂S, Jq_all, Kq_all)
end

"""
    RHFhess(ints::IntegralHelper{Float64,<:AbstractDFERI}, wfn::RHF)

Full analytic DF-RHF Hessian (Cartesian, `3*Natoms x 3*Natoms`, atomic
units). See `RHFhess.jl` for the exact-ERI dispatch.
"""
function RHFhess(ints::Fermi.Integrals.IntegralHelper{Float64,<:Fermi.Integrals.AbstractDFERI}, wfn::RHF)
    molecule = wfn.molecule
    atoms = molecule.atoms
    natm = length(atoms)

    cache = build_df_hess_cache(wfn, ints)
    bset = cache.bset
    nbas = bset.nbas

    ndocc = wfn.ndocc
    C = wfn.orbitals.C
    Co = C[:, 1:ndocc]
    Cv = C[:, ndocc+1:end]
    D = Co * Co'
    P = 2.0 .* D
    H0mat = ints["T"] .+ ints["V"]
    F = zeros(nbas, nbas)
    build_fock!(F, H0mat, D, ints)
    Q = 2.0 .* D * F * D

    responses = [_rhf_hess_response_df(wfn, ints, cache, iA, Co, Cv) for iA in 1:natm]

    H = zeros(3 * natm, 3 * natm)

    for iA in 1:natm, iB in iA:natm
        RA = responses[iA]
        RB = responses[iB]

        H2e = ERI_hess_JK_df(cache, P, iA, iB)
        H1e_S = GaussianBasis.∇2overlap(bset, iA, iB)
        H1e_T = GaussianBasis.∇2kinetic(bset, iA, iB)
        H1e_V = GaussianBasis.∇2nuclear(bset, iA, iB)
        Hnn = Molecules.∇2nuclear_repulsion(atoms, iA, iB)

        block = zeros(3, 3)
        for qA in 1:3, qB in 1:3
            dP_B = 2.0 .* (@view RB.dD[:, :, qB])
            dQ_B = 2.0 .* ((@view RB.dD[:, :, qB]) * F * D .+ D * (@view RB.dF[:, :, qB]) * D .+ D * F * (@view RB.dD[:, :, qB]))

            term = sum(dP_B .* (@view RA.∂H[:, :, qA])) - sum(dQ_B .* (@view RA.∂S[:, :, qA]))
            term += 2 * sum(dP_B .* (@view RA.Jq[:, :, qA])) - sum(dP_B .* (@view RA.Kq[:, :, qA]))

            term += sum(P .* (@view (H1e_T .+ H1e_V)[:, :, qA, qB]))
            term -= sum(Q .* (@view H1e_S[:, :, qA, qB]))
            term += H2e[qA, qB]
            term += Hnn[qA, qB]

            block[qA, qB] = term
        end

        IA = (3*(iA-1)+1):(3*iA)
        IB = (3*(iB-1)+1):(3*iB)
        H[IA, IB] .= block
        if iA != iB
            H[IB, IA] .= block'
        end
    end

    return H
end
