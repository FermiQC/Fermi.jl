# Shared first-leg transform for every compute_*! method below: builds the
# partially-contracted dense intermediate 𝐗[p,ν,ρ,σ] = Σ_μ (μν|ρσ)*C1[μ,p]
# from a SparseERI list, exploiting the same 8-fold permutation symmetry the
# original per-orbital loop did. `C1` selects which orbital set (`Co`/`Cv`)
# transforms the first leg -- occupied for OOOO/OOOV/OOVV/OVOV/OVVV, virtual
# for VVVV.
#
# This replaces an earlier version that spawned one task per column of `C1`
# and had EACH task rescan the *entire* `eri_vals` array -- i.e. O(n1 × nnz)
# redundant index-decoding and branching, not just O(n1 × nnz) FMA work
# (which is unavoidable). Measured on benzene/cc-pvdz (114 bf, 21 occ): the
# old OVOV transform took ~74s; a single pass over `eri_vals`, decoding each
# entry once and scattering into all n1 columns at once, is what `build_fock!`'s
# SparseERI method already does for the same reason (chunked worker-pool,
# private per-task accumulator, merged at the end -- see that function for
# the same pattern with a smaller, 2-index accumulator).
function _sparse_first_leg_transform(eri_vals::Vector{T}, idxs, C1::AbstractMatrix{T}, nbf::Int) where T<:AbstractFloat
    n1 = size(C1, 2)
    # Transposed so C1t[:,μ] is a contiguous length-n1 column (C1[μ,:] would
    # be a stride-nbf view) -- matters since this is read on every scatter.
    C1t = collect(C1')

    chunksize = clamp(cld(length(eri_vals), 50*Threads.nthreads()), 1, 1000)
    ntasks = Threads.nthreads()

    chunks = [collect(chunk) for chunk in Iterators.partition(eachindex(eri_vals), chunksize)]
    requests = Channel{Vector{Int}}(length(chunks))
    for c in chunks
        put!(requests, c)
    end
    close(requests)

    results = Channel{Array{T,4}}(ntasks)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            buf = zeros(T, n1, nbf, nbf, nbf)
            for chunk in requests
                @inbounds @fastmath for z in chunk
                    V = eri_vals[z]
                    μ,ν,ρ,σ = idxs[z] .+ 1
                    γμν = μ !== ν
                    γρσ = ρ !== σ
                    γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)

                    @views begin
                        Cμ = C1t[:,μ]
                        if γab && γμν && γρσ
                            Cν = C1t[:,ν]; Cρ = C1t[:,ρ]; Cσ = C1t[:,σ]
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                            buf[:,ν,σ,ρ] .+= V .* Cμ
                            buf[:,μ,ρ,σ] .+= V .* Cν
                            buf[:,μ,σ,ρ] .+= V .* Cν
                            buf[:,σ,μ,ν] .+= V .* Cρ
                            buf[:,σ,ν,μ] .+= V .* Cρ
                            buf[:,ρ,μ,ν] .+= V .* Cσ
                            buf[:,ρ,ν,μ] .+= V .* Cσ

                        elseif γab && γμν
                            Cν = C1t[:,ν]; Cρ = C1t[:,ρ]
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                            buf[:,μ,ρ,σ] .+= V .* Cν
                            buf[:,σ,μ,ν] .+= V .* Cρ
                            buf[:,σ,ν,μ] .+= V .* Cρ

                        elseif γab && γρσ
                            Cρ = C1t[:,ρ]; Cσ = C1t[:,σ]
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                            buf[:,ν,σ,ρ] .+= V .* Cμ
                            buf[:,σ,μ,ν] .+= V .* Cρ
                            buf[:,ρ,μ,ν] .+= V .* Cσ

                        elseif γμν && γρσ
                            Cν = C1t[:,ν]
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                            buf[:,ν,σ,ρ] .+= V .* Cμ
                            buf[:,μ,ρ,σ] .+= V .* Cν
                            buf[:,μ,σ,ρ] .+= V .* Cν

                        elseif γab
                            Cρ = C1t[:,ρ]
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                            buf[:,σ,μ,ν] .+= V .* Cρ

                        else
                            buf[:,ν,ρ,σ] .+= V .* Cμ
                        end
                    end
                end
            end
            put!(results, buf)
        end
    end

    out = zeros(T, n1, nbf, nbf, nbf)
    for _ in 1:ntasks
        out .+= take!(results)
    end
    return out
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    idxs = aoints["ERI"].indexes
    eri_vals = aoints["ERI"].data

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Co = T.(Co)
    end

    iνρσ = _sparse_first_leg_transform(eri_vals, idxs, Co, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        OOOO[i,j,k,l] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
    end
    I["OOOO"] = OOOO
end


function compute_OOOV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    iνρσ = _sparse_first_leg_transform(eri_vals, idxs, Co, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        OOOV[i,j,k,a] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
    end
    I["OOOV"] = OOOV
end


function compute_OOVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    iνρσ = _sparse_first_leg_transform(eri_vals, idxs, Co, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        OOVV[i,j,a,b] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
    end
    I["OOVV"] = OOVV
end


function compute_OVOV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    iνρσ = _sparse_first_leg_transform(eri_vals, idxs, Co, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        OVOV[i,a,j,b] :=  iνρσ[i, ν, ρ, σ]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    iνρσ = _sparse_first_leg_transform(eri_vals, idxs, Co, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        OVVV[i,a,b,c] := iνρσ[i, ν, ρ, σ]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
    end
    I["OVVV"] = OVVV
end


function compute_VVVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
    end

    aνρσ = _sparse_first_leg_transform(eri_vals, idxs, Cv, nbf)

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin
        VVVV[a,b,c,d] :=  aνρσ[a, ν, ρ, σ]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end
