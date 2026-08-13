# Two-electron (Coulomb + exchange) contribution to the UHF Hessian,
# contracted against FIXED densities (Dtot,Dα,Dβ) -- the "integral
# response"/direct piece only, UCPHF's density-response contribution is
# assembled separately. Direct generalization of
# RHF/Hessians/TwoElectronHess.jl to two spin densities, same relationship
# UHFgrad.jl's two-electron term has to RHFgrad.jl's -- see that file's
# header comment for the canonical-quadruple/permutation-collapse
# derivation (unchanged, depends only on the ERI's own symmetry) and
# RHF/Hessians/TwoElectronHess.jl's for the C1/C2 8-arrangement-collapse
# derivation (also unchanged for the same reason).
#
# X^{AB}_direct = 0.5*Dtot⊗Dtot*(μν|λσ)^AB - 0.5*[Dα⊗Dα+Dβ⊗Dβ]*(μν|λσ)^AB,
# matching UHFgrad.jl's X^A one derivative order up exactly the way RHF's
# ERI_hess_JK matches RHFgrad.jl's X^A.

function contract_C1_uhf(Dtot, I, J, K, L, blk)
    Dtot_ij = @view Dtot[I, J]; Dtot_kl = @view Dtot[K, L]
    c = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        c += Dtot_ij[a,b]*Dtot_kl[cc,d]*blk[a,b,cc,d]
    end
    return c
end

function contract_C2_canonical(Dα, Dβ, I, J, K, L, blk)
    Dα_il = @view Dα[I, L]; Dα_kj = @view Dα[K, J]
    Dβ_il = @view Dβ[I, L]; Dβ_kj = @view Dβ[K, J]
    ea = 0.0
    eb = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        ea += Dα_il[a,d]*Dα_kj[cc,b]*v
        eb += Dβ_il[a,d]*Dβ_kj[cc,b]*v
    end
    return ea + eb
end

function contract_C2_cross(Dα, Dβ, I, J, K, L, blk)
    Dα_jl = @view Dα[J, L]; Dα_ki = @view Dα[K, I]
    Dβ_jl = @view Dβ[J, L]; Dβ_ki = @view Dβ[K, I]
    ea = 0.0
    eb = 0.0
    @inbounds for d in axes(blk,4), cc in axes(blk,3), b in axes(blk,2), a in axes(blk,1)
        v = blk[a,b,cc,d]
        ea += Dα_jl[b,d]*Dα_ki[cc,a]*v
        eb += Dβ_jl[b,d]*Dβ_ki[cc,a]*v
    end
    return ea + eb
end

"""
    ERI_hess_JK(bset::BasisSet, Dtot, Dα, Dβ, iA::Int, iB::Int; σvals=nothing)

Two-electron Coulomb+exchange contribution to the UHF Hessian block
`(iA,iB)`, `0.5*C1 - 0.5*C2`. `Dtot,Dα,Dβ` are treated as fixed (the
converged SCF densities); the UCPHF density-response contribution is
assembled separately. `σvals` (optional, `GaussianBasis.schwarz_bounds(bset)`)
is atom-pair-independent -- callers looping over many `(iA,iB)` pairs
should compute it once and pass it through. See `RHF/Hessians/TwoElectronHess.jl`
for the RHF version this generalizes.
"""
function ERI_hess_JK(bset::BasisSet, Dtot::AbstractMatrix, Dα::AbstractMatrix, Dβ::AbstractMatrix, iA::Int, iB::Int; σvals = nothing)
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

                        C1_local[x, y] += mult * contract_C1_uhf(Dtot, I, J, K, L, blk)

                        if i == j && k == l
                            Vcan_C2 = contract_C2_canonical(Dα, Dβ, I, J, K, L, blk)
                            C2_local[x, y] += ij_ne_kl ? 2*Vcan_C2 : Vcan_C2
                        else
                            Vcan_C2 = contract_C2_canonical(Dα, Dβ, I, J, K, L, blk)
                            Vcross_C2 = contract_C2_cross(Dα, Dβ, I, J, K, L, blk)
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

    return 0.5 .* C1 .- 0.5 .* C2
end
