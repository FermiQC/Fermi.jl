# CPHF (Coupled-Perturbed Hartree-Fock) orbital-response solver.
#
# The occ-virt rotation U^(y) (y = a nuclear Cartesian perturbation) that
# describes how the MOs respond to a geometry displacement solves
#
#   (eps_a - eps_i) U_ai + sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj = -B_ai
#
# The coupling-matrix action A*U (the sum over b,j) is the only genuinely new
# two-electron contraction this needs, and it turns out to be exactly a
# standard Fock-build contraction in disguise -- no new ERI-handling code,
# just a different (symmetrized, trial) density plugged into the existing
# build_fock! dispatches.
#
# Derivation: let ΔD_μν = sum_ai [C_μa U_ai C_νi + C_μi U_ai C_νa] (the
# symmetric AO-basis density built from the trial rotation U, shape
# (nvir,ndocc)), and ΔF_μν = 2*ΔD_rs*(μν|rs) - ΔD_rs*(μr|νs) -- exactly
# build_fock!'s two-electron contraction, with ΔD standing in for the density
# and no one-electron piece added. Expanding ΔD in terms of U and transforming
# ΔF to the MO occ-virt block:
#
#   sum_μν C_μa ΔF_μν C_νi
#     = 4*sum_bj U_bj*(ai|bj) - sum_bj U_bj*[(ab|ij) + (aj|ib)]
#     = sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj
#
# using (μν|jb)=(μν|bj) (AO ERI symmetry under the ket-pair swap) to combine
# the two Coulomb-shaped terms, and the same relabeling for the exchange
# piece -- exactly the CPHF coupling matrix's action on U. So:
# Cv'*ΔF*Co = A*U, with ΔF obtained by literally reusing build_fock!.

"""
    build_response_fock!(ΔF, U, Co, Cv, ints)

Two-electron "response Fock" contraction used as the CPHF orbital-Hessian
matrix-vector product A*U. `U` is a trial occupied-virtual rotation, shape
`(nvir,ndocc)`; `Co`, `Cv` are the occupied/virtual blocks of the MO
coefficient matrix. See this file's header comment for the derivation.
"""
function build_response_fock!(ΔF, U, Co, Cv, ints)
    ΔD = Cv * U * Co'
    ΔD .+= ΔD'
    H0 = zeros(size(ΔF))
    build_fock!(ΔF, H0, ΔD, ints)
    return ΔF
end

# Scratch-buffer core: ΔD/H0/tmp are caller-owned, reused across calls
# instead of allocating fresh. This is the CPHF matrix-vector product,
# called once per CG iteration (up to cphf_max_iter times) per Cartesian
# direction per atom, so the old allocating form's fresh zeros(nbas,nbas)
# added up the same way the gradient's/direct-Hessian's per-quartet
# allocation did (see RHFgrad.jl's and GaussianBasis.jl's scratch-buffer
# docstrings) -- smaller in absolute magnitude here, same avoidable churn.
function build_response_fock!(ΔF, U, Co, Cv, ints, ΔD, H0, tmp)
    mul!(tmp, Cv, U)
    mul!(ΔD, tmp, Co')
    ΔD .+= ΔD'
    build_fock!(ΔF, H0, ΔD, ints)
    return ΔF
end

"""
    cphf_Amatvec(U, Co, Cv, ints)

The CPHF coupling matrix's action on a trial occ-virt rotation `U`,
`sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj`, returned as a `(nvir,ndocc)`
matrix in the same shape/index convention as `U`. Thin wrapper around
`build_response_fock!` for use as a KrylovKit linear map.
"""
function cphf_Amatvec(U, Co, Cv, ints)
    nbas = size(Co, 1)
    ΔF = zeros(nbas, nbas)
    build_response_fock!(ΔF, U, Co, Cv, ints)
    return Cv' * ΔF * Co
end

# Scratch-buffer core of cphf_Amatvec, called once per CG iteration --
# result/ΔF/tmp2 are caller-owned alongside build_response_fock!'s own
# ΔD/H0/tmp, same motivation as that function's scratch-buffer overload.
function cphf_Amatvec!(result, U, Co, Cv, ints, ΔF, ΔD, H0, tmp, tmp2)
    build_response_fock!(ΔF, U, Co, Cv, ints, ΔD, H0, tmp)
    mul!(tmp2, Cv', ΔF)
    mul!(result, tmp2, Co)
    return result
end

# --- CPHF right-hand side ---
#
# The CPHF equation is (εa-εi)U_ai + (A·U)_ai = -B_ai. B has two pieces:
#
#   1. The "skeleton" Fock derivative: differentiate the AO integrals
#      directly (H, ERI) with the density D=Co*Co' held FIXED at its
#      converged value -- exactly the quantity RHFgrad.jl already builds one
#      contraction short of (it goes on to dot this against P,Q to get a
#      gradient; here it's kept as a full AO matrix and transformed to the
#      MO occ-virt block instead), minus ε_i*S_ai^(y) from the generalized
#      eigenvalue problem's metric.
#
#   2. A correction from how the occupied-occupied block of the MO
#      coefficient response is fixed (not solved for) by the orthonormality
#      constraint C(y)'S(y)C(y)=I: writing C_i^(y) = sum_p C_p U_pi, the
#      symmetric part of U is forced to U_pq+U_qp=-S_pq^(y) (MO-transformed
#      overlap derivative); the standard gauge choice U_ij=-0.5*S_ij^(y) for
#      occ-occ i,j fixes that block outright (a pure occ-occ rotation
#      wouldn't change the density; this symmetric piece isn't a rotation,
#      so it DOES contribute to dD/dy through a known, computable term,
#      dD_S/dy = -Co*S_oo^(y)*Co', S_oo^(y):=Co'S^(y)Co). Since F depends on
#      the density D exactly the same way build_response_fock! already
#      exploits, this piece's contribution to B is just
#      build_fock!(0, dD_S/dy, ints) transformed to the ai block -- same
#      reuse trick, no new two-electron code.

"""
    eri_grad_JK(bset::BasisSet, D::AbstractMatrix, iA::Int; ij_vals=nothing, σvals=nothing)

Coulomb- and exchange-shaped contractions of the *first*-derivative ERI
against a density-like matrix `D`, for all three Cartesian directions of
atom `iA`: `Jq[m,n] = sum_rs D[r,s]*d(mn|rs)/dA_q`, `Kq[m,n] = sum_rs
D[r,s]*d(mr|ns)/dA_q`. Returned as two `(nbas,nbas,3)` arrays.

`ij_vals`/`σvals` (`GaussianBasis.schwarz_bounds(bset)`) are optional and
atom-independent -- callers looping over many atoms should compute them
once and pass them through rather than recomputing per call.
"""
function eri_grad_JK(bset::BasisSet, D::AbstractMatrix, iA::Int; ij_vals = nothing, σvals = nothing)
    # Built directly from GaussianBasis.∇sparseERI_2e4c's compressed
    # (permutation-unique) derivative list -- never materializes the dense
    # (nbas,nbas,nbas,nbas,3) tensor ∇ERI_2e4c! would. For each stored
    # unique quartet, accumulates into Jq/Kq over the (deduplicated)
    # 8-member index-permutation orbit that shares its value; validated
    # against the dense ∇ERI_2e4c!-based reference to machine precision
    # (sto-3g and cc-pVDZ, several atoms) before use here.
    nbas = bset.nbas
    Jq = zeros(nbas, nbas, 3)
    Kq = zeros(nbas, nbas, 3)
    idx, xyz... = GaussianBasis.∇sparseERI_2e4c(bset, iA; ij_vals=ij_vals, σvals=σvals)
    for i in eachindex(idx)
        μ, ν, λ, σ = Int.(idx[i])
        orbit = unique([(μ, ν, λ, σ), (ν, μ, λ, σ), (μ, ν, σ, λ), (ν, μ, σ, λ),
                        (λ, σ, μ, ν), (σ, λ, μ, ν), (λ, σ, ν, μ), (σ, λ, ν, μ)])
        for q in 1:3
            ∇k = xyz[q][i]
            abs(∇k) < 1e-12 && continue
            for (a, b, c, d) in orbit
                Jq[a, b, q] += D[c, d] * ∇k
                Kq[a, c, q] += D[b, d] * ∇k
            end
        end
    end
    return Jq, Kq
end

"""
    cphf_rhs(wfn, ints, iA)

CPHF right-hand side `B_ai^(y)` for all three Cartesian directions of atom
`iA`, returned as a `(nvir,ndocc,3)` array, alongside the `∂H,∂S,Jq_all,
Kq_all` intermediates it was built from (`(nbas,nbas,3)` each).
"""
function cphf_rhs(wfn::RHF, ints, iA; ij_vals = nothing, σvals = nothing)
    # ∂H/∂S/Jq_all/Kq_all are returned (not just B) because
    # _rhf_hess_response needs those same intermediates again for its own
    # Fskel/dF assembly -- returning them here instead of forcing a second,
    # identical ∇kinetic!/∇nuclear!/∇overlap!/eri_grad_JK pass saved ~97s
    # of ~276s total CPHF time on a real 24-atom molecule (see
    # cphf_solve_full's comment).
    bset = BasisSet(wfn.orbitals.basis, wfn.molecule.atoms)
    nbas = bset.nbas
    ndocc = wfn.ndocc
    nvir = nbas - ndocc
    C = wfn.orbitals.C
    Co = @view C[:, 1:ndocc]
    Cv = @view C[:, ndocc+1:end]
    eps_o = wfn.orbitals.eps[1:ndocc]
    D = Co * Co'

    ∂H = zeros(nbas, nbas, 3)
    ∂Vtmp = zeros(nbas, nbas, 3)
    GaussianBasis.∇kinetic!(∂H, bset, iA)
    GaussianBasis.∇nuclear!(∂Vtmp, bset, iA)
    ∂H .+= ∂Vtmp

    ∂S = zeros(nbas, nbas, 3)
    GaussianBasis.∇overlap!(∂S, bset, iA)

    Jq_all, Kq_all = eri_grad_JK(bset, D, iA; ij_vals=ij_vals, σvals=σvals)

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
    cphf_solve(wfn, ints, iA)

Solve the CPHF equations for the occupied-virtual orbital response
`U_ai^(y)` to displacing atom `iA`, for all three Cartesian directions,
returned as a `(nvir,ndocc,3)` array. Thin wrapper around `cphf_solve_full`
for callers that only need `U` -- see that function for the solve itself.
"""
function cphf_solve(wfn::RHF, ints, iA; ij_vals = nothing, σvals = nothing)
    U, _, _, _, _ = cphf_solve_full(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)
    return U
end

"""
    cphf_solve_full(wfn, ints, iA)

Same solve as `cphf_solve`, but also returns `cphf_rhs`'s intermediate
`∂H,∂S,Jq_all,Kq_all` alongside `U`, for callers (like `_rhf_hess_response`)
that need those too. Solves `(εa-εi)U + cphf_Amatvec(U) = -cphf_rhs(...)`
via `KrylovKit.linsolve` with (Jacobi-preconditioned) CG.
"""
function cphf_solve_full(wfn::RHF, ints, iA; ij_vals = nothing, σvals = nothing)
    # Returning ∂H/∂S/Jq_all/Kq_all here (not just U) lets _rhf_hess_response
    # reuse them instead of recomputing via a second ∇kinetic!/∇nuclear!/
    # ∇overlap!/eri_grad_JK pass -- profiling a real 24-atom molecule found
    # that redundant pass costing ~97s of ~276s total CPHF time.
    #
    # Preconditioning: plain CG on (εa-εi).*x + cphf_Amatvec(x) was found
    # (via profiling) to hit cphf_max_iter (default 50) on ~85% of
    # right-hand sides for a real molecule -- an operator-conditioning
    # problem (the εa-εi denominators span a huge range, from near-
    # degenerate valence pairs to deep-core/high-virtual gaps of tens of
    # Hartree), not an algorithm bug. Since M := Δeps.*x + cphf_Amatvec(x)
    # is SPD and Δeps is diagonal (trivially invertible), the standard fix
    # (used by essentially every CPHF/TDDFT iterative solver, including
    # PySCF's Z-vector CPHF) is to solve the congruent, still-SPD system
    # D^-1/2 M D^-1/2 y = D^-1/2 b (D:=Diagonal(vec(Δeps))) instead --
    # substituting x = y./d turns this into
    # y .+ cphf_Amatvec(y./d)./d = -B./d (elementwise ./d, no need to
    # materialize D or its inverse), recovering U = y./d at the end.
    ndocc = wfn.ndocc
    nbas = size(wfn.orbitals.C, 1)
    nvir = nbas - ndocc
    C = wfn.orbitals.C
    Co = @view C[:, 1:ndocc]
    Cv = @view C[:, ndocc+1:end]
    eps = wfn.orbitals.eps
    Δeps = eps[ndocc+1:end] .- eps[1:ndocc]'  # (nvir,ndocc), εa - εi
    d = sqrt.(Δeps)                           # Jacobi/diagonal preconditioner scaling

    maxiter = Options.get("cphf_max_iter")
    tol = Options.get("cphf_conv")

    B, ∂H, ∂S, Jq_all, Kq_all = cphf_rhs(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)
    U = zeros(nvir, ndocc, 3)

    # Reused across every CG iteration (all 3 directions) instead of
    # letting cphf_Amatvec/build_response_fock! allocate ΔF/ΔD/H0/tmp/tmp2
    # fresh on every matvec call -- see cphf_Amatvec!'s docstring.
    ΔF = zeros(nbas, nbas)
    ΔD = zeros(nbas, nbas)
    H0 = zeros(nbas, nbas)
    tmp = zeros(nbas, ndocc)
    tmp2 = zeros(nvir, nbas)
    result = zeros(nvir, ndocc)

    for q in 1:3
        function precond_matvec(y)
            cphf_Amatvec!(result, y ./ d, Co, Cv, ints, ΔF, ΔD, H0, tmp, tmp2)
            return y .+ result ./ d
        end
        rhs = -B[:, :, q] ./ d
        sol, info = KrylovKit.linsolve(precond_matvec, rhs, zeros(nvir, ndocc), KrylovKit.CG(; maxiter=maxiter, tol=tol))
        U[:, :, q] .= sol ./ d
    end
    return U, ∂H, ∂S, Jq_all, Kq_all
end
