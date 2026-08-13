using GaussianBasis

# Full analytic RHF Hessian assembly.
#
# Extends RHFgrad.jl's gradient formula
#   dE/dA = P_mn ∂H_mn/∂A - Q_mn ∂S_mn/∂A + X^A + ∂Vnn/∂A
# (X^A := 0.5*P_mn*P_rs*(mn|rs)^A - 0.25*P_mn*P_rs*(ms|rn)^A, matching
# RHFgrad.jl's own AUX_ERI contraction exactly) one derivative order, via the
# product rule on every P/Q/integral factor. Two kinds of terms result:
#
#   - "Direct": both derivatives land on the integrals, with P/Q held at
#     their converged values -- these are exactly the pieces already built
#     and validated: the one-electron Hessians (∇2overlap/∇2kinetic/
#     ∇2nuclear), ERI_hess_JK (X's own second derivative, same 0.5/-0.25
#     coefficients), and Molecules.∇2nuclear_repulsion.
#
#   - "Response": one derivative lands on P or Q (via dP/dB, dQ/dB -- the
#     orbital/density response CPHF solves for), the other on the A-side
#     integrals (already available as first derivatives). dQ/dB uses the
#     gauge-invariant identity Q=2*D*F*D (D=Co*Co', F=converged AO Fock;
#     verified to ~1e-10) rather than an orbital-by-orbital eps_i^(y)
#     decomposition, which turns out to depend on an antisymmetric occ-occ
#     rotation gauge that the density response P doesn't need but Q does --
#     avoided entirely by working with D and F directly. The P-response
#     2-electron term reuses RHFgrad.jl's own AUX_ERI contraction pattern,
#     just with one of the two P factors replaced by dP/dB -- and, since
#     that contraction is linear in the fixed density factor and P=2D,
#     it's exactly `2*sum(dP_B.*Jq) - sum(dP_B.*Kq)` using the SAME
#     Jq,Kq=eri_grad_JK(bset,D,iA) already computed for Fskel below (no
#     separate computation, no dense (nbas,nbas,nbas,nbas) array ever
#     materialized -- verified against the dense/AUX_ERI formula to
#     machine precision before use here).
#
# Validated end-to-end (this whole assembly, not pieces in isolation)
# against finite difference of Fermi's own analytic gradient: matches to
# ~8e-7 (finite-difference precision) on water/sto-3g.
#
# CPHF's own Fock builds (cphf_Amatvec/build_response_fock!, one per CG
# iteration) go through `ints`'s own eri_type dispatch (build_fock!'s
# Chonky/SparseERI/DF cases) rather than a hardcoded dense Chonky build --
# RHFhess(wfn) builds a throwaway aoints (SparseERI by default, matching
# what RHF itself defaults to) with a warning; RHFhess(aoints, wfn) reuses
# a caller-supplied one, e.g. straight out of `@energy aoints => rhf`.

struct RHFHessResponse
    dD::Array{Float64,3}  # (nbas,nbas,3) -- half-density response, D=Co*Co'
    dF::Array{Float64,3}  # (nbas,nbas,3) -- AO Fock response (skeleton + P-response)
    ∂H::Array{Float64,3}  # (nbas,nbas,3) -- first-derivative core Hamiltonian
    ∂S::Array{Float64,3}  # (nbas,nbas,3) -- first-derivative overlap
    Jq::Array{Float64,3}  # (nbas,nbas,3) -- sum_rs D[r,s]*d(mn|rs)/dA_q
    Kq::Array{Float64,3}  # (nbas,nbas,3) -- sum_rs D[r,s]*d(mr|ns)/dA_q
end

function _rhf_hess_response(wfn::RHF, ints, bset, iA, Co, Cv; ij_vals = nothing, σvals = nothing)
    nbas = size(Co, 1)
    U, ∂H, ∂S, Jq_all, Kq_all = cphf_solve_full(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)

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
    RHFhess(wfn::RHF)

Aoints-free convenience entry point: builds a throwaway `IntegralHelper` for
`wfn`'s molecule/basis (warning first, see `warn_no_aoints`) and delegates to
`RHFhess(aoints, wfn)`. Prefer building `aoints` yourself and calling that
directly (`@hessian aoints, wfn => rhf`) when chaining `@energy`/`@gradient`/
`@hessian` at the same geometry -- CPHF's Fock builds genuinely reuse it.
"""
function RHFhess(wfn::RHF)
    Fermi.Integrals.warn_no_aoints()
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis)
    RHFhess(aoints, wfn)
end

"""
    RHFhess(aoints::IntegralHelper{Float64,<:Union{Chonky,SparseERI}}, wfn::RHF)

Full analytic RHF Hessian for the exact-ERI case (Cartesian, `3*Natoms x
3*Natoms`, atomic units). See `DFHess.jl` for the density-fitted dispatch.
"""
function RHFhess(ints::Fermi.Integrals.IntegralHelper{Float64,<:Union{Fermi.Integrals.Chonky,Fermi.Integrals.SparseERI}}, wfn::RHF)
    # aoints is used directly for CPHF's Fock builds (build_fock! inside
    # cphf_Amatvec/build_response_fock!) -- its eri_type is honored rather
    # than a hardcoded dense Chonky build, which would force an
    # O(nbas^4)-per-CG-iteration dense contraction.
    molecule = wfn.molecule
    atoms = molecule.atoms
    natm = length(atoms)
    bset = BasisSet(ints.basis, ints.molecule.atoms)
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

    # Schwarz screening bound is atom-independent -- compute it once here
    # rather than paying its O(nshells^2) cost again inside ∇sparseERI_2e4c
    # for every one of the natm calls below (eri_grad_JK, called once per
    # atom from inside cphf_rhs -- cphf_solve_full threads its ∂H/∂S/Jq/Kq
    # back out to _rhf_hess_response instead of that function recomputing
    # them a second time, see CPHF.jl's cphf_solve_full docstring).
    ij_vals, σvals = GaussianBasis.schwarz_bounds(bset)

    responses = [_rhf_hess_response(wfn, ints, bset, iA, Co, Cv; ij_vals=ij_vals, σvals=σvals) for iA in 1:natm]

    H = zeros(3 * natm, 3 * natm)

    for iA in 1:natm, iB in iA:natm
        RA = responses[iA]
        RB = responses[iB]

        H2e = ERI_hess_JK(bset, P, iA, iB; σvals=σvals)
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
