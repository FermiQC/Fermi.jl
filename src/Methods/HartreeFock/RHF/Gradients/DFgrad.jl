using GaussianBasis
using Molecules
using TensorOperations
using LinearAlgebra

# Analytic RHF gradient for a density-fitted (JKFIT/RIFIT) wave function.
#
# The two-electron energy contribution X^A = 0.5*P_mn*P_rs*(mn|rs) -
# 0.25*P_mn*P_rs*(ms|rn) (same X^A as RHFgrad.jl's Chonky/Sparse dispatches)
# is rewritten with the RI approximation (mn|rs) ~ sum_PQ (mn|P) J^-1_PQ
# (Q|rs), J_PQ := (P|Q) the 2-center auxiliary metric. No CPHF/orbital
# response is needed here, same as the exact-ERI gradient case -- the SCF
# procedure is variational in D regardless of whether the two-electron
# integrals it was converged against are exact or RI-approximated, so
# Hellmann-Feynman still applies: D (and the auxiliary basis, which simply
# moves rigidly with its assigned atom like any other basis function) can be
# held fixed and only the integrals differentiated.
#
# Coulomb (J) part: with c_P := P_mn*(mn|P) and d_P := J^-1_PQ*c_Q,
#   0.5*P_mn*P_rs*(mn|rs) = 0.5*c.d
#   d/dA = dc.d - 0.5*d.(dJ/dA).d      (using d(J^-1)/dA = -J^-1 (dJ/dA) J^-1)
#
# Exchange (K) part needs an occupied-orbital half-transform to stay at DF's
# usual cost (exactly the trick DirectHelper.jl's exchange_to_fock! already
# uses for the energy): W[m,i,P] := sum_n C_ni*(mn|P), V[m,i,P] :=
# sum_Q W[m,i,Q]*J^-1_QP. Differentiating the exchange energy expression
# (E_K := -sum_mn,i,PQ D_mn*W[m,i,P]*J^-1_PQ*W[n,i,Q], matching
# exchange_to_fock!'s own -1 sign convention) and using
# P_ms*P_rn*(ms|rn) = -4*E_K gives, for the X^A contribution
# (-0.25 * d[-4*E_K]/dA = dE_K/dA):
#   dE_K/dA = -2*sum_mn,i D_mn*dW[m,i,P]*J^-1_PQ*W[n,i,Q]  +  Tr(Qtot.dJ/dA)
# where Qtot_PQ := sum_i (D*V_i)'*V_i (a small Naux x Naux matrix,
# accumulated one occupied orbital at a time to avoid ever materializing an
# (nbas,nbas,Naux,Naux) intermediate).

"""
    RHFgrad(aoints::IntegralHelper{Float64,<:AbstractDFERI}, wfn::RHF)

Analytic RHF gradient for a density-fitted (`JKFIT`/`RIFIT`) wave function --
see this file's header for the RI-HF gradient formula. Full 3-center/2-center
derivative tensors are materialized per atom (mirroring `RHFgrad`'s Chonky
dispatch), which is tractable here since DF tensors scale as `nbas^2*naux`,
not `nbas^4`.
"""
function RHFgrad(aoints::Fermi.Integrals.IntegralHelper{Float64,<:Fermi.Integrals.AbstractDFERI}, wfn::RHF)

    atoms = wfn.molecule.atoms
    Natoms = length(atoms)
    bset = BasisSet(aoints.basis, aoints.molecule.atoms)
    auxbset = aoints.eri_type.basisset
    nbas = bset.nbas
    naux = auxbset.nbas
    ndocc = wfn.ndocc

    ∂E = zeros(Natoms, 3)

    @views Co = wfn.orbitals.C[:, 1:ndocc]
    D = Co * Co'
    P = 2.0 * D
    Q = 2.0 * Co * diagm(wfn.orbitals.eps[1:ndocc]) * Co'

    Pmn = GaussianBasis.ERI_2e3c(bset, auxbset)   # (nbas,nbas,naux)
    Jmat = GaussianBasis.ERI_2e2c(auxbset)        # (naux,naux)
    Jinv = inv(Jmat)

    @tensor c[A] := P[m, n] * Pmn[m, n, A]
    d = Jinv * c

    @tensor W[m, i, A] := Co[n, i] * Pmn[m, n, A]   # (nbas,ndocc,naux)
    @tensor V[m, i, A] := W[m, i, B] * Jinv[B, A]   # (nbas,ndocc,naux)

    # Qtot = sum_i V_i' * D * V_i, accumulated one occupied orbital at a time
    Qtot = zeros(naux, naux)
    for i in 1:ndocc
        Vi = @view V[:, i, :]
        Qtot .+= (D * Vi)' * Vi
    end

    ∂H = zeros(nbas, nbas, 3)
    ∂S = zeros(nbas, nbas, 3)
    ∂Vtmp = zeros(nbas, nbas, 3)
    ∇Pmn = zeros(nbas, nbas, naux, 3)
    ∇Jmat = zeros(naux, naux, 3)
    AUX = zeros(nbas, nbas)

    for a in eachindex(atoms)

        ∂H .= 0
        ∂S .= 0
        ∂Vtmp .= 0
        ∇Pmn .= 0
        ∇Jmat .= 0

        ∂E[a, :] .= Molecules.∇nuclear_repulsion(atoms, a)

        GaussianBasis.∇kinetic!(∂H, bset, a)
        GaussianBasis.∇nuclear!(∂Vtmp, bset, a)
        ∂H .+= ∂Vtmp

        GaussianBasis.∇overlap!(∂S, bset, a)

        GaussianBasis.∇ERI_2e3c!(∇Pmn, bset, auxbset, a)
        GaussianBasis.∇ERI_2e2c!(∇Jmat, auxbset, a)

        for q in 1:3
            @views vH = ∂H[:, :, q]
            AUX .= P .* vH
            ∂E[a, q] += sum(AUX)

            @views vS = ∂S[:, :, q]
            AUX .= Q .* vS
            ∂E[a, q] -= sum(AUX)

            @views ∇Pq = ∇Pmn[:, :, :, q]
            @views ∇Jq = ∇Jmat[:, :, q]

            @tensor dc[A] := P[m, n] * ∇Pq[m, n, A]
            Jpart = dot(dc, d) - 0.5 * dot(d, ∇Jq * d)

            @tensor dW[m, i, A] := Co[n, i] * ∇Pq[m, n, A]
            @tensor dWJ[m, i, B] := dW[m, i, A] * Jinv[A, B]
            @tensor term1 = D[m, n] * dWJ[m, i, B] * W[n, i, B]
            Kpart = -2.0 * term1 + dot(Qtot, ∇Jq)

            ∂E[a, q] += Jpart + Kpart
        end
    end

    return ∂E
end
