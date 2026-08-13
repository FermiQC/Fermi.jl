using GaussianBasis
using Molecules

function UHFgrad(x...)
    UHFgrad(Molecule(), x...)
end

function UHFgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic" || dtype == "auto"
        # SCF itself needs atomic integrals regardless -- warn once here
        # (not inside UHF(mol,Alg), which would also build one) and reuse
        # the same aoints for the gradient call below (though the gradient
        # itself doesn't need it, see UHFgrad(wfn::UHF)'s docstring).
        Fermi.Integrals.warn_no_aoints()
        aoints = Fermi.Integrals.IntegralHelper{Float64}(molecule=mol)
        wfn = UHF(aoints)
        # x... (e.g. an explicit eri_type) falls back to the older dispatch.
        isempty(x) ? UHFgrad(aoints, wfn) : UHFgrad(wfn, x...)
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.HartreeFock.UHF, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

"""
    UHFgrad(wfn::UHF)

Analytic UHF gradient. Builds a throwaway `IntegralHelper` for `wfn`'s
molecule/basis and delegates to the `(aoints, wfn)` dispatch matching its
`eri_type` -- mirrors `RHFgrad(wfn::RHF)` (`RHF/Gradients/RHFgrad.jl`)
exactly, see that docstring for why no `warn_no_aoints` call is needed here.
"""
function UHFgrad(wfn::UHF)
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis)
    UHFgrad(aoints, wfn)
end

"""
    UHFgrad(wfn::UHF, eri_type::Fermi.Integrals.AbstractERI)

Explicit-`eri_type` entry point, mirrors `RHFgrad(wfn::RHF,eri_type)`'s
docstring/rationale exactly (avoids the same infinite-recursion trap that
motivated adding that method).
"""
function UHFgrad(wfn::UHF, eri_type::Fermi.Integrals.AbstractERI)
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis, eri_type=eri_type)
    UHFgrad(aoints, wfn)
end

# --- Two-electron contribution: integral-direct, shell-quartet-by-shell-quartet ---
#
# Direct generalization of RHFgrad.jl's two-electron gradient to two spin
# densities -- see that file's header comment for the full integral-direct/
# canonical-quadruple/permutation-collapse derivation, which applies
# unchanged here (it depends only on the ERI's own permutation symmetry,
# not on which density gets contracted against it).
#
# UHF energy (from UHFHelper.jl's UHFEnergy/build_fock!, Dtot:=Dα+Dβ,
# Jtot:=Jα+Jβ, Fα=H+Jtot-Kα, Fβ=H+Jtot-Kβ): E = Dtot·H + 0.5*Jtot·Dtot -
# 0.5*(Kα·Dα+Kβ·Dβ) + Vnn. UHF is variational in Dα,Dβ independently, so the
# same Hellmann-Feynman argument RHFgrad.jl's X^A uses gives (holding
# Dα,Dβ fixed, differentiating only the explicit integral dependence):
#
#   X^A = 0.5*Dtot_μν*Dtot_λσ*(μν|λσ)^A - 0.5*[Dα_μλ*Dα_νσ + Dβ_μλ*Dβ_νσ]*(μν|λσ)^A
#
# Verified this reduces EXACTLY to RHF's X^A = 0.5*P·P·(mn|rs)^A -
# 0.25*P·P·(ms|rn)^A in the closed-shell limit (Dα=Dβ=D, P=2D, Dtot=2D):
# 0.5*4*D·D*(mn|rs)^A - 0.5*2*D·D*(ms|rn)^A = 2*D·D*(mn|rs)^A - D·D*(ms|rn)^A,
# matching RHF's 0.5*4*D·D*(mn|rs)^A - 0.25*4*D·D*(ms|rn)^A term-for-term --
# this is exploited directly in test_UHF.jl as a free, strong correctness
# check (UHF forced closed-shell must match the already-validated RHF
# gradient to near machine precision, independent of any finite-difference
# or Psi4 comparison).
#
# `contract_canonical`/`contract_cross`/`quartet_contribution` below share
# their names with RHFgrad.jl's own (different arity -- three densities
# instead of one -- so Julia's normal multiple dispatch disambiguates,
# exactly like UHFHelper.jl's `build_fock!` already coexists with RHF's).

function contract_canonical(Dtot, Dα, Dβ, I, J, K, L, blk)
    Dtot_ij = @view Dtot[I, J]; Dtot_kl = @view Dtot[K, L]
    Dα_il = @view Dα[I, L]; Dα_kj = @view Dα[K, J]
    Dβ_il = @view Dβ[I, L]; Dβ_kj = @view Dβ[K, J]
    c = 0.0
    ea = 0.0
    eb = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        c += Dtot_ij[a,b]*Dtot_kl[cc,d]*v
        ea += Dα_il[a,d]*Dα_kj[cc,b]*v
        eb += Dβ_il[a,d]*Dβ_kj[cc,b]*v
    end
    return 0.5*c - 0.5*(ea+eb)
end

function contract_cross(Dtot, Dα, Dβ, I, J, K, L, blk)
    Dtot_ij = @view Dtot[I, J]; Dtot_kl = @view Dtot[K, L]
    Dα_jl = @view Dα[J, L]; Dα_ki = @view Dα[K, I]
    Dβ_jl = @view Dβ[J, L]; Dβ_ki = @view Dβ[K, I]
    c = 0.0
    ea = 0.0
    eb = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        c += Dtot_ij[a,b]*Dtot_kl[cc,d]*v
        ea += Dα_jl[b,d]*Dα_ki[cc,a]*v
        eb += Dβ_jl[b,d]*Dβ_ki[cc,a]*v
    end
    return 0.5*c - 0.5*(ea+eb)
end

function quartet_contribution(Dtot, Dα, Dβ, I, J, K, L, blk_q, i, j, k, l)
    Vcan = contract_canonical(Dtot, Dα, Dβ, I, J, K, L, blk_q)
    ij_ne_kl = GaussianBasis.index2(i-1,j-1) != GaussianBasis.index2(k-1,l-1)
    if i == j && k == l
        return ij_ne_kl ? 2*Vcan : Vcan
    else
        Vcross = contract_cross(Dtot, Dα, Dβ, I, J, K, L, blk_q)
        n_indep = (i != j ? 2 : 1) * (k != l ? 2 : 1) * (ij_ne_kl ? 2 : 1)
        half = n_indep ÷ 2
        return half*Vcan + half*Vcross
    end
end

"""
    UHFgrad(aoints::IntegralHelper{Float64,<:Union{Chonky,SparseERI}}, wfn::UHF)

Analytic UHF gradient for the exact-ERI case. Returns the Cartesian
gradient as a `(Natoms, 3)` matrix, atomic units. See `RHFgrad.jl` for the
RHF version this generalizes.
"""
function UHFgrad(aoints::Fermi.Integrals.IntegralHelper{Float64,<:Union{Fermi.Integrals.Chonky,Fermi.Integrals.SparseERI}}, wfn::UHF)

    atoms = wfn.molecule.atoms
    Natoms = length(atoms)
    bset = BasisSet(aoints.basis, aoints.molecule.atoms)
    nbas = bset.nbas
    nshells = bset.nshells

    ∂E = zeros(Natoms, 3)

    Nα = wfn.molecule.Nα
    Nβ = wfn.molecule.Nβ
    @views Coα = wfn.orbitals.Cα[:,1:Nα]
    @views Coβ = wfn.orbitals.Cβ[:,1:Nβ]

    Dα = Coα * Coα'
    Dβ = Coβ * Coβ'
    Dtot = Dα .+ Dβ

    Q = Coα * diagm(wfn.orbitals.epsα[1:Nα]) * Coα' .+ Coβ * diagm(wfn.orbitals.epsβ[1:Nβ]) * Coβ'

    # Preallocate arrays
    ∂H = zeros(nbas, nbas, 3)
    ∂S = zeros(nbas, nbas, 3)

    # Intermediate auxiliary arrays
    AUX = zeros(nbas, nbas)

    Nvals = GaussianBasis.num_basis.(bset.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:nshells]
    Nmax = maximum(Nvals)

    # Schwarz screening bound is atom-independent -- compute it once here
    # rather than paying its O(nshells^2) cost again for every atom below.
    _, σvals = GaussianBasis.schwarz_bounds(bset)
    cutoff = 1e-12

    for a in eachindex(atoms)

        ∂H .= 0
        ∂S .= 0

        # Nuclear repulsion
        ∂E[a, :] .= Molecules.∇nuclear_repulsion(atoms, a)

        # Store kinetic into H and nuclear into S
        GaussianBasis.∇kinetic!(∂H, bset, a)
        GaussianBasis.∇nuclear!(∂S, bset, a)

        ∂H .+= ∂S
        ∂S .= 0

        # Now use S for overlap
        GaussianBasis.∇overlap!(∂S, bset, a)

        for q in 1:3
            @views vH = ∂H[:,:,q]
            AUX .= Dtot .* vH
            ∂E[a, q] += sum(AUX)

            @views vS = ∂S[:,:,q]
            AUX .= Q .* vS
            ∂E[a, q] -= sum(AUX)
        end

        A = atoms[a]
        on_atom = falses(nshells)
        for s in 1:nshells
            bset.basis[s].atom === A && (on_atom[s] = true)
        end

        # Pre-screen the canonical (i≤j,k≤l,ij-pair≥kl-pair) quartet work
        # list -- Schwarz bound and atom-membership are both O(1) checks,
        # cheap to do up front so worker tasks only ever see quartets that
        # need real work.
        quartets = NTuple{4,Int}[]
        for i in 1:nshells, j in i:nshells
            σij = σvals[GaussianBasis.index2(i-1,j-1)+1]
            ij_idx = GaussianBasis.index2(i-1,j-1)
            for k in 1:nshells, l in k:nshells
                ij_idx < GaussianBasis.index2(k-1,l-1) && continue
                σkl = σvals[GaussianBasis.index2(k-1,l-1)+1]
                σij*σkl <= cutoff && continue

                i_A, j_A, k_A, l_A = on_atom[i], on_atom[j], on_atom[k], on_atom[l]
                (!(i_A||j_A||k_A||l_A) || (i_A&&j_A&&k_A&&l_A)) && continue

                push!(quartets, (i, j, k, l))
            end
        end

        isempty(quartets) && continue

        ntasks = Threads.nthreads()
        chunksize = clamp(cld(length(quartets), 4*ntasks), 1, 200)
        requests = Channel{Vector{NTuple{4,Int}}}(Inf)
        for chunk in Iterators.partition(quartets, chunksize)
            put!(requests, collect(chunk))
        end
        close(requests)

        results = Channel{Vector{Float64}}(ntasks)

        @sync for _ in 1:ntasks
            Threads.@spawn begin
                local_E = zeros(3)
                # Per-task scratch, reused across every quartet this worker
                # ever sees -- see RHFgrad.jl's identical comment/rationale.
                out_buf = zeros(Nmax, Nmax, Nmax, Nmax, 3)
                cint_buf = zeros(Cdouble, 3*Nmax^4)
                cint_tmp = zeros(Cdouble, 3*Nmax^4)
                shls_buf = zeros(Cint, 4)
                for chunk in requests
                    for (i, j, k, l) in chunk
                        Ni, Nj, Nk, Nl = Nvals[i], Nvals[j], Nvals[k], Nvals[l]
                        I = (ao_offset[i]+1):(ao_offset[i]+Ni)
                        J = (ao_offset[j]+1):(ao_offset[j]+Nj)
                        K = (ao_offset[k]+1):(ao_offset[k]+Nk)
                        L = (ao_offset[l]+1):(ao_offset[l]+Nl)

                        on_A = (on_atom[i], on_atom[j], on_atom[k], on_atom[l])
                        blk = @view out_buf[1:Ni, 1:Nj, 1:Nk, 1:Nl, :]
                        GaussianBasis.∇ERI_2e4c!(blk, bset, on_A, i, j, k, l, cint_buf, cint_tmp, shls_buf)

                        for q in 1:3
                            bq = @view blk[:,:,:,:,q]
                            local_E[q] += quartet_contribution(Dtot, Dα, Dβ, I, J, K, L, bq, i, j, k, l)
                        end
                    end
                end
                put!(results, local_E)
            end
        end

        for _ in 1:ntasks
            ∂E[a, :] .+= take!(results)
        end
    end

    return ∂E
end
