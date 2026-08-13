# Two-electron (Coulomb + exchange) contribution to the RHF Hessian, contracted
# against a FIXED density matrix P -- i.e. the "integral response"/direct piece
# only (mirrors RHFgrad.jl's existing gradient formula one derivative order up;
# CPHF's density-*response* contribution, needing dP/dR, is a separate term
# assembled elsewhere). Never materializes the raw (nbas,nbas,nbas,nbas,3,3) ERI
# Hessian tensor -- contracts against P immediately, shell-quartet by
# shell-quartet, same spirit as a Fock build.
#
# GaussianBasis.jl's ∇2ERI_2e4c(bset,iA,iB,i,j,k,l) is the light-bridge
# shell-quartet Hessian primitive this contracts against -- the raw libcint
# kernel plumbing this file used to do directly (cint2e_ipip1/ipvip1/ip1ip2,
# with their non-obvious kernel-to-shell-pair mapping and reversed axis
# orientation for cross-shell kernels) now lives there instead, ported
# verbatim rather than touched, since getting that right the first time
# already needed careful cint_funcs.h reading and finite-difference checks.
#
# Visits only CANONICAL shell quadruples (i≤j, k≤l, (i,j)-pair ≥ (k,l)-pair --
# the same restriction RHFgrad.jl's outer loop uses), one `∇2ERI_2e4c` call
# each. Getting BOTH C1 and C2's full contribution (summed over all 8
# permutation-symmetric arrangements of that quadruple) from that SINGLE call
# needed its own derivation, verified against a brute-force reference (each
# of the 8 arrangements computed independently via its own permuted view of
# the canonical block, exactly as `RHFgrad.jl`'s own `quartet_contribution`
# derivation was verified) over 1000 random trials before trusting it:
#
#   C1 (Coulomb-shape, weight P[I,J]⊗P[K,L]) is IDENTICAL across all 8
#   arrangements -- same result RHFgrad.jl's Coulomb term found. So a
#   canonical quadruple's total C1 contribution is just `mult * Vcan_C1`,
#   `mult` being the count of distinct arrangements (1, 2, 4, or 8 depending
#   on i==j?/k==l?/(i,j)-pair==(k,l)-pair? degeneracies -- same combinatorics
#   as RHFgrad.jl's `quartet_contribution`).
#
#   C2 (exchange-shape, weight P[I,L]⊗P[K,J]) splits into exactly the same
#   "canonical-class"/"cross-class" (4-and-4, or fewer when degenerate)
#   pattern RHFgrad.jl's own exchange term did -- `half*Vcan_C2 +
#   half*Vcross_C2` in general, `Vcan_C2` alone (with the right multiplier)
#   when i==j and k==l simultaneously (bra-swap and ket-swap both become
#   no-ops, so there's no "cross" arrangement to distinguish).
#
# A pleasant side effect: since C1 and C2 now both derive from the SAME
# canonical block (i,j,k,l) -- not two independent calls, a genuinely
# different integral needing its own Schwarz bound, as an earlier version of
# this file had to handle -- one Schwarz screening decision (σ_ij·σ_kl) now
# correctly covers both terms.
#
# Three efficiency fixes found by profiling a real (24-atom caffeine, not a
# toy 3-atom molecule) system, none of which are visible at toy scale:
#
#   1. `σvals` is now an optional kwarg, computed once by RHFhess.jl's main
#      loop and threaded through, instead of every one of the O(natm^2)
#      ERI_hess_JK calls recomputing the same O(nshells^2) Schwarz bound
#      pass from scratch (RHFgrad.jl already avoided this; this file hadn't
#      been updated to match).
#
#   2. Atom membership is checked via a precomputed BitVector (`on_A`/`on_B`,
#      O(nshells), built once per call) instead of `any(p -> ...
#      shellval[p]].atom == Aat, 1:4)` -- a closure re-evaluated on every one
#      of the ~nshells^4/8 canonical quadruples. Profiling a single
#      ERI_hess_JK call on caffeine/sto-3g (52 shells, so ~950k canonical
#      quadruples to filter) showed this closure-based check alone consuming
#      ~26% of total samples -- more than either of the two real libcint
#      kernels it was gating. RHFgrad.jl's loop already used the cheaper
#      precomputed-array form; this file hadn't been brought in line with it.
#
#   3. Threaded via the same channel-based worker-pool pattern RHFgrad.jl
#      uses (and, originally, RHFHelper.jl's SparseERI build_fock!): work
#      items (surviving canonical quadruples) go into a bounded Channel,
#      each Threads.@spawn'd worker drains chunks into its own local (3,3)
#      accumulators, summed on the main task only after every worker
#      finishes. This file had no threading at all until now.

function contract_C1(P, I, J, K, L, blk)
    Pij = @view P[I, J]; Pkl = @view P[K, L]
    c = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        c += Pij[a,b]*Pkl[cc,d]*blk[a,b,cc,d]
    end
    return c
end

function contract_C2_canonical(P, I, J, K, L, blk)
    Pil = @view P[I, L]; Pkj = @view P[K, J]
    e = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        e += Pil[a,d]*Pkj[cc,b]*blk[a,b,cc,d]
    end
    return e
end

function contract_C2_cross(P, I, J, K, L, blk)
    Pjl = @view P[J, L]; Pki = @view P[K, I]
    e = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        e += Pjl[b,d]*Pki[cc,a]*blk[a,b,cc,d]
    end
    return e
end

"""
    ERI_hess_JK(bset::BasisSet, P::AbstractMatrix, iA::Int, iB::Int; σvals=nothing)

Two-electron Coulomb+exchange contribution to the RHF Hessian block
`(iA,iB)`, `0.5*C1 - 0.25*C2` where `C1[x,y] = sum_mnrs
P[m,n]*P[r,s]*d2(mn|rs)/dA_x dB_y` (Coulomb-shape) and `C2[x,y] = sum_mnrs
P[m,n]*P[r,s]*d2(ms|rn)/dA_x dB_y` (exchange-shape). `P` is treated as
fixed (the converged SCF density); the CPHF density-response contribution
is assembled separately. `σvals` (optional, `GaussianBasis.schwarz_bounds(bset)`)
is atom-pair-independent -- callers looping over many `(iA,iB)` pairs
should compute it once and pass it through.
"""
function ERI_hess_JK(bset::BasisSet, P::AbstractMatrix, iA::Int, iB::Int; σvals = nothing)
    Nvals = GaussianBasis.num_basis.(bset.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:bset.nshells]
    nshells = bset.nshells
    Nmax = maximum(Nvals)

    if σvals === nothing
        _, σvals = GaussianBasis.schwarz_bounds(bset)
    end
    cutoff = 1e-12

    Aat = bset.atoms[iA]
    Bat = bset.atoms[iB]
    on_A = falses(nshells)
    on_B = falses(nshells)
    for s in 1:nshells
        bset.basis[s].atom === Aat && (on_A[s] = true)
        bset.basis[s].atom === Bat && (on_B[s] = true)
    end

    # Pre-screen the canonical quadruple work list -- Schwarz bound and
    # atom-membership are both O(1) checks (given precomputed on_A/on_B),
    # cheap to do up front so worker tasks only ever see quadruples that
    # need real work. Also excludes the iA==iB, all-four-shells-on-that-atom
    # translational-invariance case here (GaussianBasis.jl's ∇2ERI_2e4c no
    # longer checks it itself when called with precomputed flags below --
    # this loop already has everything needed to check it once, up front).
    quartets = NTuple{4,Int}[]
    @inbounds for i in 1:nshells, j in i:nshells
        ij_idx = GaussianBasis.index2(i-1,j-1)
        σij = σvals[ij_idx+1]
        for k in 1:nshells, l in k:nshells
            ij_idx < GaussianBasis.index2(k-1,l-1) && continue
            σkl = σvals[GaussianBasis.index2(k-1,l-1)+1]
            σij*σkl <= cutoff && continue

            (on_A[i]||on_A[j]||on_A[k]||on_A[l]) || continue
            (on_B[i]||on_B[j]||on_B[k]||on_B[l]) || continue
            (iA == iB && on_A[i] && on_A[j] && on_A[k] && on_A[l]) && continue

            push!(quartets, (i, j, k, l))
        end
    end

    C1 = zeros(3, 3)
    C2 = zeros(3, 3)

    isempty(quartets) && return C1

    ntasks = Threads.nthreads()
    chunksize = clamp(cld(length(quartets), 4*ntasks), 1, 200)
    requests = Channel{Vector{NTuple{4,Int}}}(Inf)
    for chunk in Iterators.partition(quartets, chunksize)
        put!(requests, collect(chunk))
    end
    close(requests)

    results = Channel{Tuple{Matrix{Float64},Matrix{Float64}}}(ntasks)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            C1_local = zeros(3, 3)
            C2_local = zeros(3, 3)
            # Per-task scratch, reused across every quartet this worker ever
            # sees -- mirrors RHFgrad.jl's fix for the same allocation
            # pattern (see GaussianBasis.jl's scratch-buffer ∇2ERI_2e4c!
            # docstring); this loop calls it up to 16x per quartet (4x4
            # posA/posB combinations) versus the gradient's 4 branches, so
            # the old allocating path cost proportionally more here.
            out_buf = zeros(Nmax, Nmax, Nmax, Nmax, 3, 3)
            cint_buf = zeros(Cdouble, 9*Nmax^4)
            cint_t1 = zeros(Cdouble, 9*Nmax^4)
            cint_t2 = zeros(Cdouble, 9*Nmax^4)
            shls_buf = zeros(Cint, 4)
            for chunk in requests
                for (i, j, k, l) in chunk
                    Ni, Nj, Nk, Nl = Nvals[i], Nvals[j], Nvals[k], Nvals[l]
                    ioff, joff, koff, loff = ao_offset[i], ao_offset[j], ao_offset[k], ao_offset[l]
                    I = (ioff+1):(ioff+Ni)
                    J = (joff+1):(joff+Nj)
                    K = (koff+1):(koff+Nk)
                    L = (loff+1):(loff+Nl)

                    Xflag = (on_A[i], on_A[j], on_A[k], on_A[l])
                    Yflag = (on_B[i], on_B[j], on_B[k], on_B[l])
                    blk_full = @view out_buf[1:Ni, 1:Nj, 1:Nk, 1:Nl, :, :]
                    GaussianBasis.∇2ERI_2e4c!(blk_full, bset, Xflag, Yflag, i, j, k, l, cint_buf, cint_t1, cint_t2, shls_buf)

                    ij_ne_kl = GaussianBasis.index2(i-1,j-1) != GaussianBasis.index2(k-1,l-1)
                    mult = (i != j ? 2 : 1) * (k != l ? 2 : 1) * (ij_ne_kl ? 2 : 1)

                    for y in 1:3, x in 1:3
                        blk = @view blk_full[:,:,:,:,x,y]

                        C1_local[x, y] += mult * contract_C1(P, I, J, K, L, blk)

                        if i == j && k == l
                            Vcan_C2 = contract_C2_canonical(P, I, J, K, L, blk)
                            C2_local[x, y] += ij_ne_kl ? 2*Vcan_C2 : Vcan_C2
                        else
                            Vcan_C2 = contract_C2_canonical(P, I, J, K, L, blk)
                            Vcross_C2 = contract_C2_cross(P, I, J, K, L, blk)
                            half = mult ÷ 2
                            C2_local[x, y] += half*Vcan_C2 + half*Vcross_C2
                        end
                    end
                end
            end
            put!(results, (C1_local, C2_local))
        end
    end

    for _ in 1:ntasks
        C1_local, C2_local = take!(results)
        C1 .+= C1_local
        C2 .+= C2_local
    end

    return 0.5 .* C1 .- 0.25 .* C2
end
