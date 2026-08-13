using GaussianBasis

# Full analytic UHF Hessian assembly. Direct generalization of
# RHF/Hessians/RHFhess.jl to two spin densities -- same "direct" (both
# derivatives on integrals, densities fixed) / "response" (one derivative
# via UCPHF's dD/dB, other on first-derivative A-side integrals) split, via
# the product rule on UHFgrad.jl's gradient formula
#   dE/dA = Dtot·∂H_A - Qtot·∂S_A + X^A(A) + ∂Vnn_A
# one derivative order further, w.r.t. a second atom B:
#
#   d/dB[Dtot·∂H_A] = dDtot_B·∂H_A + Dtot·∂²H_AB
#   d/dB[Qtot·∂S_A] = dQtot_B·∂S_A + Qtot·∂²S_AB
#   d/dB[X^A(A)]    = [dDtot_B·Jq_tot(A) - dDα_B·Kqα(A) - dDβ_B·Kqβ(A)]  (response)
#                     + ERI_hess_JK(A,B)                                  (direct)
#
# (the response piece's coefficients are all 1, not RHF's 2 -- Dtot has no
# RHF-style P=2D doubling; same reason UHFgrad.jl's X^A has 0.5/0.5 instead
# of RHF's 0.5/0.25).
#
# Qα (energy-weighted density response) uses the same gauge-invariant
# identity RHF's Q=2*D*F*D does, one spin at a time and without the factor
# of 2: at self-consistency Fα*Cα=S*Cα*diag(εα) gives Dα*Fα*Dα =
# Coα*(Coα'*Fα*Coα)*Coα' = Coα*(Coα'*S*Coα*diag(εα_o))*Coα' =
# Coα*diag(εα_o)*Coα' = Qα exactly (Coα'*S*Coα=I, Coα orthonormal). So
# dQα_B = dDα_B*Fα*Dα + Dα*dFα_B*Dα + Dα*Fα*dDα_B, same three-term product
# rule RHF's dQ_B uses, per spin.

struct UHFHessResponse
    dDα::Array{Float64,3}  # (nbas,nbas,3) -- density response, spin α
    dDβ::Array{Float64,3}  # (nbas,nbas,3) -- density response, spin β
    dFα::Array{Float64,3}  # (nbas,nbas,3) -- AO Fock response, spin α (skeleton + response)
    dFβ::Array{Float64,3}  # (nbas,nbas,3) -- AO Fock response, spin β
    ∂H::Array{Float64,3}   # (nbas,nbas,3) -- first-derivative core Hamiltonian
    ∂S::Array{Float64,3}   # (nbas,nbas,3) -- first-derivative overlap
    Jq::Array{Float64,3}   # (nbas,nbas,3) -- sum_rs Dtot[r,s]*d(mn|rs)/dA_q
    Kqα::Array{Float64,3}  # (nbas,nbas,3) -- sum_rs Dα[r,s]*d(mr|ns)/dA_q
    Kqβ::Array{Float64,3}  # (nbas,nbas,3) -- sum_rs Dβ[r,s]*d(mr|ns)/dA_q
end

function _uhf_hess_response(wfn::UHF, ints, bset, iA, Coα, Cvα, Coβ, Cvβ; ij_vals = nothing, σvals = nothing)
    nbas = size(Coα, 1)
    (Uα, Uβ), ∂H, ∂S, Jq_tot, Kqα, Kqβ = ucphf_solve_full(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)

    dDα = zeros(nbas, nbas, 3)
    dDβ = zeros(nbas, nbas, 3)
    dFα = zeros(nbas, nbas, 3)
    dFβ = zeros(nbas, nbas, 3)
    for q in 1:3
        # No factor of 2 on Jq_tot -- see UCPHF.jl's identical fix/note.
        Fskelα = @view(∂H[:, :, q]) .+ (@view Jq_tot[:, :, q]) .- (@view Kqα[:, :, q])
        Fskelβ = @view(∂H[:, :, q]) .+ (@view Jq_tot[:, :, q]) .- (@view Kqβ[:, :, q])

        Sq = @view ∂S[:, :, q]
        S_oo_α = Coα' * Sq * Coα
        S_oo_β = Coβ' * Sq * Coβ
        dD_S_α = -Coα * S_oo_α * Coα'
        dD_S_β = -Coβ * S_oo_β * Coβ'
        Gα = zeros(nbas, nbas)
        Gβ = zeros(nbas, nbas)
        _response_fock_from_densities!(Gα, Gβ, dD_S_α, dD_S_β, ints)

        ΔFα = zeros(nbas, nbas)
        ΔFβ = zeros(nbas, nbas)
        build_response_fock!(ΔFα, ΔFβ, Uα[:, :, q], Uβ[:, :, q], Coα, Cvα, Coβ, Cvβ, ints)

        dDqα = Cvα * Uα[:, :, q] * Coα'
        dDqα .+= dDqα'
        dDqα .+= dD_S_α
        dDα[:, :, q] .= dDqα

        dDqβ = Cvβ * Uβ[:, :, q] * Coβ'
        dDqβ .+= dDqβ'
        dDqβ .+= dD_S_β
        dDβ[:, :, q] .= dDqβ

        dFα[:, :, q] .= Fskelα .+ ΔFα .+ Gα
        dFβ[:, :, q] .= Fskelβ .+ ΔFβ .+ Gβ
    end

    return UHFHessResponse(dDα, dDβ, dFα, dFβ, ∂H, ∂S, Jq_tot, Kqα, Kqβ)
end

"""
    UHFhess(wfn::UHF)

Aoints-free convenience entry point, mirrors `RHFhess(wfn::RHF)` exactly.
"""
function UHFhess(wfn::UHF)
    Fermi.Integrals.warn_no_aoints()
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis)
    UHFhess(aoints, wfn)
end

"""
    UHFhess(aoints::IntegralHelper{Float64,<:Union{Chonky,SparseERI}}, wfn::UHF)

Full analytic UHF Hessian (Cartesian, `3*Natoms x 3*Natoms`, atomic units).
See `RHFhess.jl` for the RHF version this generalizes, and this file's
header comment for the term-by-term derivation.
"""
function UHFhess(ints::Fermi.Integrals.IntegralHelper{Float64,<:Union{Fermi.Integrals.Chonky,Fermi.Integrals.SparseERI}}, wfn::UHF)
    molecule = wfn.molecule
    atoms = molecule.atoms
    natm = length(atoms)
    bset = BasisSet(ints.basis, ints.molecule.atoms)
    nbas = bset.nbas

    Nα = molecule.Nα
    Nβ = molecule.Nβ
    Cα = wfn.orbitals.Cα
    Cβ = wfn.orbitals.Cβ
    Coα = Cα[:, 1:Nα]; Cvα = Cα[:, Nα+1:end]
    Coβ = Cβ[:, 1:Nβ]; Cvβ = Cβ[:, Nβ+1:end]
    Dα = Coα * Coα'
    Dβ = Coβ * Coβ'
    Dtot = Dα .+ Dβ

    Jα = zeros(nbas, nbas); Jβ = zeros(nbas, nbas)
    Kα = zeros(nbas, nbas); Kβ = zeros(nbas, nbas)
    Fα = zeros(nbas, nbas); Fβ = zeros(nbas, nbas)
    build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints)

    Qα = Dα * Fα * Dα
    Qβ = Dβ * Fβ * Dβ
    Qtot = Qα .+ Qβ

    # Schwarz screening bound is atom-independent -- compute it once here,
    # see RHFhess.jl's identical comment/rationale (eri_grad_JK_uhf, called
    # once per atom from inside ucphf_rhs -- ucphf_solve_full threads its
    # ∂H/∂S/Jq/Kqα/Kqβ back out to _uhf_hess_response instead of
    # recomputing them a second time).
    ij_vals, σvals = GaussianBasis.schwarz_bounds(bset)

    responses = [_uhf_hess_response(wfn, ints, bset, iA, Coα, Cvα, Coβ, Cvβ; ij_vals=ij_vals, σvals=σvals) for iA in 1:natm]

    H = zeros(3 * natm, 3 * natm)

    for iA in 1:natm, iB in iA:natm
        RA = responses[iA]
        RB = responses[iB]

        H2e = ERI_hess_JK(bset, Dtot, Dα, Dβ, iA, iB; σvals=σvals)
        H1e_S = GaussianBasis.∇2overlap(bset, iA, iB)
        H1e_T = GaussianBasis.∇2kinetic(bset, iA, iB)
        H1e_V = GaussianBasis.∇2nuclear(bset, iA, iB)
        Hnn = Molecules.∇2nuclear_repulsion(atoms, iA, iB)

        block = zeros(3, 3)
        for qA in 1:3, qB in 1:3
            dDα_B = @view RB.dDα[:, :, qB]
            dDβ_B = @view RB.dDβ[:, :, qB]
            dDtot_B = dDα_B .+ dDβ_B

            dFα_B = @view RB.dFα[:, :, qB]
            dFβ_B = @view RB.dFβ[:, :, qB]

            dQα_B = dDα_B * Fα * Dα .+ Dα * dFα_B * Dα .+ Dα * Fα * dDα_B
            dQβ_B = dDβ_B * Fβ * Dβ .+ Dβ * dFβ_B * Dβ .+ Dβ * Fβ * dDβ_B
            dQtot_B = dQα_B .+ dQβ_B

            term = sum(dDtot_B .* (@view RA.∂H[:, :, qA])) - sum(dQtot_B .* (@view RA.∂S[:, :, qA]))
            term += sum(dDtot_B .* (@view RA.Jq[:, :, qA])) - sum(dDα_B .* (@view RA.Kqα[:, :, qA])) - sum(dDβ_B .* (@view RA.Kqβ[:, :, qA]))

            term += sum(Dtot .* (@view (H1e_T .+ H1e_V)[:, :, qA, qB]))
            term -= sum(Qtot .* (@view H1e_S[:, :, qA, qB]))
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
