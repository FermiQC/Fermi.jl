using GaussianBasis
using Molecules

function RHFgrad(x...)
    RHFgrad(Molecule(), x...)
end

function RHFgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic" || dtype == "auto"
        # SCF itself needs atomic integrals regardless -- warn once here
        # (not inside RHF(mol,Alg), which would also build one) and reuse
        # the same aoints for the gradient call below (though the gradient
        # itself doesn't need it, see RHFgrad(wfn::RHF)'s docstring).
        Fermi.Integrals.warn_no_aoints()
        aoints = Fermi.Integrals.IntegralHelper{Float64}(molecule=mol)
        wfn = RHF(aoints)
        # x... (e.g. an explicit eri_type) falls back to the older dispatch.
        isempty(x) ? RHFgrad(aoints, wfn) : RHFgrad(wfn, x...)
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.HartreeFock.RHF, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

"""
    RHFgrad(wfn::RHF)

Analytic RHF gradient. Builds a throwaway `IntegralHelper` for `wfn`'s
molecule/basis and delegates to the dispatch matching its `eri_type`
(exact-ERI, this file; density-fitted, `DFgrad.jl`), following the current
`@set df <bool>` option.
"""
function RHFgrad(wfn::RHF)
    # No warn_no_aoints call here, on purpose: unlike RHFhess, no dispatch
    # reached from here touches aoints's J/K cache, only cheap basis/
    # molecule metadata -- so there's no redundant integral computation to
    # warn about. Matters for cost too, since this runs every geometry-
    # optimizer iteration and IntegralHelper construction itself is cheap.
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis)
    RHFgrad(aoints, wfn)
end

"""
    RHFgrad(wfn::RHF, eri_type::Fermi.Integrals.AbstractERI)

Analytic RHF gradient with an explicit `eri_type` (e.g. `@gradient rhf <=
Fermi.Integrals.Chonky()`), overriding the current `@set df` option.
"""
function RHFgrad(wfn::RHF, eri_type::Fermi.Integrals.AbstractERI)
    # Needed to avoid an infinite loop: without this method, RHFgrad(wfn,x...)
    # falls through to the fully generic RHFgrad(x...) catch-all instead of
    # dispatching here, which re-enters RHFgrad(mol::Molecule,x...),
    # reconverges a FRESH SCF, and recurses indefinitely. (Reached via
    # @gradient's macro-expanded string template, not a literal call site --
    # easy to miss with a plain grep, which is how this method got
    # accidentally deleted once before.)
    aoints = Fermi.Integrals.IntegralHelper(molecule=wfn.molecule, basis=wfn.orbitals.basis, eri_type=eri_type)
    RHFgrad(aoints, wfn)
end

# --- Two-electron contribution: integral-direct, shell-quartet-by-shell-quartet ---
#
# Unlike energy integrals (reused across every SCF/CPHF iteration, so
# materializing/caching them pays for itself), a derivative integral is only
# ever needed ONCE -- to help form X^A -- and then discarded. There is no
# reuse being given up by never storing one, so "compute a quartet, contract
# it immediately, discard it" has no downside here, unlike integral-direct
# SCF's usual compute-vs-memory tradeoff. This replaces both an older
# dense-array method (materialized the full (nbas,nbas,nbas,nbas,3) tensor
# via ∇ERI_2e4c!, now deleted -- dead code, nothing called it) and, until
# now, a "sparse array" method that still fully built a compressed
# (idx,∇x,∇y,∇z) list for the whole atom via ∇sparseERI_2e4c before
# contracting it in a second pass. Neither materialization step is needed:
# GaussianBasis.jl's shell-quartet-level ∇ERI_2e4c(bset,iA,i,j,k,l) primitive
# is exactly the building block for going one step further.
#
# X^A = 0.5*P_μν*P_λσ*(μν|λσ)^A - 0.25*P_μν*P_λσ*(μσ|λν)^A, summed
# unrestricted over all μ,ν,λ,σ. Decomposing that unrestricted sum by shell
# quadruple, a SINGLE shell quadruple (I,J,K,L)'s contribution (using
# blk[a,b,c,d] := ∂(μ_aν_b|λ_cσ_d)/∂R, μ_a∈I etc.) is
#
#   0.5*Σ P[I,J].*P[K,L].*blk - 0.25*Σ P[I,L].*P[K,J].*blk        ("canonical")
#
# (`contract_canonical` below) -- valid for ANY shell assignment (I,J,K,L),
# no ordering assumed. Visiting every (i,j,k,l) shell quadruple unrestricted
# would need up to 8x more `∇ERI_2e4c` calls than necessary; this instead
# loops only over CANONICAL unique quadruples -- i≤j, k≤l, AND the (i,j)
# shell-pair index ≥ the (k,l) one -- the same restriction `∇ERI_2e4c!`'s
# own whole-array loop uses to build its `unique_idx` list.
#
# What's still needed from a single canonical quadruple is the SUM of the
# (up to 8) symmetric arrangements' own contributions. A first attempt
# computed each of the 7 "mirror" arrangements explicitly via `permutedims`
# on the already-computed `blk` (no new libcint calls, reusing the identity
# that raw integral values are unconditionally symmetric under bra-swap
# ((μν|λσ)=(νμ|λσ)), ket-swap ((μν|λσ)=(μν|σλ)), and pair-swap
# ((μν|λσ)=(λσ|μν)) -- validated against a fresh Psi4 reference,
# water/cc-pvtz, to ~7e-10). That worked, but profiling showed
# TensorOperations/Strided's per-call dispatch overhead (StridedView setup,
# backend selection, `permutedims`'s own allocation) completely dominating
# the runtime for the very small (often 1x1x1x1-ish) blocks typical here --
# the 8 mirror computations are individually cheap but the machinery
# built for large, few tensor contractions has fixed costs that don't
# amortize at this scale.
#
# Algebraic simplification (verified against the permutedims-based version
# above over 500 random trials, max diff ~5.7e-14, before trusting it):
# working out each of the 8 arrangements' OWN P-weighted contraction in
# terms of the SAME un-permuted `blk` (substituting the permuted-index
# identity directly into the sum and relabeling dummy indices, rather than
# materializing a transposed array) shows they collapse into only TWO
# distinct values -- "canonical" (P[I,J]⊗P[K,L] Coulomb, P[I,L]⊗P[K,J]
# exchange) and "cross" (SAME Coulomb, but P[J,L]⊗P[K,I] exchange) --
# with every arrangement in {canonical, bra+ket-swapped, pair-swapped,
# pair+bra+ket-swapped} equal to "canonical", and every arrangement in
# {bra-swapped alone, ket-swapped alone, pair+bra-swapped, pair+ket-swapped}
# equal to "cross". So the full sum is just `half*Vcan + half*Vcross`,
# `half` being the count of independent axes among {i≠j, k≠l, ij-pair≠
# kl-pair} that are actually true, divided by 2 -- except when i==j AND
# k==l simultaneously, where "cross" isn't a distinct arrangement at all
# (bra-swap and ket-swap are both no-ops) and only "canonical" survives.
# `contract_canonical`/`contract_cross` use plain nested loops (not
# TensorOperations) since these blocks are too small for its machinery to
# pay for itself.
#
# Threaded via the same channel-based worker-pool pattern already used and
# validated elsewhere in this codebase (RHFHelper.jl's SparseERI
# build_fock!): work items go into a bounded Channel, each Threads.@spawn'd
# worker drains chunks into its OWN local accumulator (never indexed by
# Threads.threadid() -- exactly the pattern that was previously found unsafe
# under Julia 1.12's interactive thread pool / task migration, see
# CLAUDE.md's "Current State" notes), and workers' local accumulators are
# summed on the main task only after every worker has finished. GaussianBasis
# `∇ERI_2e4c`'s own libcint calls are already known to be safe under
# concurrent use -- `∇ERI_2e4c!`/`∇sparseERI_2e4c` already call them from
# multiple threads via GaussianBasis's own `workerpool`.

function contract_canonical(P, I, J, K, L, blk)
    Pij = @view P[I, J]; Pkl = @view P[K, L]
    Pil = @view P[I, L]; Pkj = @view P[K, J]
    c = 0.0
    e = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        c += Pij[a,b]*Pkl[cc,d]*v
        e += Pil[a,d]*Pkj[cc,b]*v
    end
    return 0.5*c - 0.25*e
end

function contract_cross(P, I, J, K, L, blk)
    Pij = @view P[I, J]; Pkl = @view P[K, L]
    Pjl = @view P[J, L]; Pki = @view P[K, I]
    c = 0.0
    e = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        c += Pij[a,b]*Pkl[cc,d]*v
        e += Pjl[b,d]*Pki[cc,a]*v
    end
    return 0.5*c - 0.25*e
end

# Sums a canonical quartet's contribution AND every other symmetric
# arrangement's, without any extra libcint calls or array permutation --
# see this file's header comment for the derivation.
function quartet_contribution(P, I, J, K, L, blk_q, i, j, k, l)
    Vcan = contract_canonical(P, I, J, K, L, blk_q)
    ij_ne_kl = GaussianBasis.index2(i-1,j-1) != GaussianBasis.index2(k-1,l-1)
    if i == j && k == l
        return ij_ne_kl ? 2*Vcan : Vcan
    else
        Vcross = contract_cross(P, I, J, K, L, blk_q)
        n_indep = (i != j ? 2 : 1) * (k != l ? 2 : 1) * (ij_ne_kl ? 2 : 1)
        half = n_indep ÷ 2
        return half*Vcan + half*Vcross
    end
end

"""
    RHFgrad(aoints::IntegralHelper{Float64,<:Union{Chonky,SparseERI}}, wfn::RHF)

Analytic RHF gradient for the exact-ERI case (`Chonky`/`SparseERI`
`eri_type`s). Returns the Cartesian gradient as a `(Natoms, 3)` matrix,
atomic units. See `DFgrad.jl` for the density-fitted dispatch, and this
file's header comment for the integral-direct implementation strategy.
"""
function RHFgrad(aoints::Fermi.Integrals.IntegralHelper{Float64,<:Union{Fermi.Integrals.Chonky,Fermi.Integrals.SparseERI}}, wfn::RHF)

    # Following eq. on C.3. Szabo & Ostlund
    atoms = wfn.molecule.atoms
    Natoms = length(atoms)
    bset = BasisSet(aoints.basis, aoints.molecule.atoms)
    nbas = bset.nbas
    nshells = bset.nshells

    ∂E = zeros(Natoms, 3)

    @views Co = wfn.orbitals.C[:,1:wfn.ndocc]

    P = 2.0 * Co * Co'
    Q = 2.0 * Co * diagm(wfn.orbitals.eps[1:wfn.ndocc]) * Co'

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
            AUX .= P .* vH
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
                # ever sees -- avoids the ~1.3 KB/call (fresh `out`, fresh
                # libcint `buf`, allocating `permutedims`) that profiling
                # found here, which added up to several GB of GC churn over
                # a full gradient. See GaussianBasis.jl's scratch-accepting
                # `∇ERI_2e4c!` docstring for the corresponding fix on its side.
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
                            local_E[q] += quartet_contribution(P, I, J, K, L, bq, i, j, k, l)
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
