# UCPHF (Unrestricted Coupled-Perturbed Hartree-Fock) orbital-response
# solver -- generalizes RHF/Hessians/CPHF.jl to two spin channels.
#
# Two coupled occ-virt rotations, Uα^(y) (nvirα,noccα) and Uβ^(y)
# (nvirβ,noccβ), solve
#
#   (εα_a-εα_i)Uα_ai + Σ_bj A^{αα}_{ai,bj}Uα_bj + Σ_bj A^{αβ}_{ai,bj}Uβ_bj = -Bα_ai
#   (εβ_a-εβ_i)Uβ_ai + Σ_bj A^{βα}_{ai,bj}Uα_bj + Σ_bj A^{ββ}_{ai,bj}Uβ_bj = -Bβ_ai
#
# with A^{αα}_{ai,bj}=(ai|bj)-(ab|ij), A^{ββ}_{ai,bj}=(ai|bj)-(ab|ij) (same
# form as RHF's coupling matrix, one spin at a time), and A^{αβ}_{ai,bj} =
# A^{βα}_{ai,bj} = (ai|bj) -- Coulomb-only cross-spin coupling, no exchange
# (UHF's Kα only ever contracts Dα, Kβ only Dβ). This is *exactly* what
# UHF's ordinary Fock build already does to any density pair (Jtot=Jα+Jβ
# shared, Kα/Kβ same-spin-only) -- so the coupled matvec is UHF's existing
# J/K machinery applied to a trial (ΔDα,ΔDβ) pair instead of the real
# density, same reuse trick RHF's CPHF.jl exploits with `build_fock!`, just
# one layer lower here since UHF's `build_fock!` (UHFHelper.jl) always adds
# the one-electron H internally (no H0-override parameter the way RHF's
# does) -- `_calcJK_uhf!` below calls UHFHelper.jl's `calcJ!`/`calcK!`/
# `calcJK!` directly (the J/K-only building blocks `build_fock!` itself
# calls) to get the same "H0=0" effect.

function _calcJK_uhf!(Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    calcJ!(Jα, Dα, ints); calcJ!(Jβ, Dβ, ints)
    calcK!(Kα, Dα, ints); calcK!(Kβ, Dβ, ints)
end

function _calcJK_uhf!(Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    calcJ!(Jα, Dα, ints); calcJ!(Jβ, Dβ, ints)
    calcK!(Kα, Dα, ints); calcK!(Kβ, Dβ, ints)
end

function _calcJK_uhf!(Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,<:SparseERI,AtomicOrbitals})
    calcJK!(Jα, Kα, Dα, ints)
    calcJK!(Jβ, Kβ, Dβ, ints)
end

# Bundled into a struct (rather than threaded as separate positional args,
# CPHF.jl's style) since UHF's doubled (α,β) channels would otherwise mean
# ~18 positional scratch arguments per call.
"""
    UCPHFScratch(nbas, nvirα, noccα, nvirβ, noccβ)

Scratch buffers for `ucphf_Amatvec!`'s call chain, reused across every CG
iteration of one `ucphf_solve_full` call instead of allocating fresh.
"""
struct UCPHFScratch
    ΔFα::Matrix{Float64}; ΔFβ::Matrix{Float64}
    ΔDα::Matrix{Float64}; ΔDβ::Matrix{Float64}
    Jα::Matrix{Float64}; Jβ::Matrix{Float64}
    Kα::Matrix{Float64}; Kβ::Matrix{Float64}
    tmpα::Matrix{Float64}; tmpβ::Matrix{Float64}     # Cvσ*Uσ, (nbas,noccσ)
    tmp2α::Matrix{Float64}; tmp2β::Matrix{Float64}   # Cvσ'*ΔFσ, (nvirσ,nbas)
end

function UCPHFScratch(nbas::Int, nvirα::Int, noccα::Int, nvirβ::Int, noccβ::Int)
    UCPHFScratch(
        zeros(nbas, nbas), zeros(nbas, nbas),
        zeros(nbas, nbas), zeros(nbas, nbas),
        zeros(nbas, nbas), zeros(nbas, nbas),
        zeros(nbas, nbas), zeros(nbas, nbas),
        zeros(nbas, noccα), zeros(nbas, noccβ),
        zeros(nvirα, nbas), zeros(nvirβ, nbas),
    )
end

"""
    _response_fock_from_densities!(ΔFα, ΔFβ, ΔDα, ΔDβ, ints)

Core J/K-only (no one-electron piece) response-Fock build shared by
`build_response_fock!` (below, trial rotation -> trial density -> this) and
`ucphf_rhs` (occ-occ orthonormality correction density -> this directly, no
rotation involved). `ΔFα = Jtot(ΔDα,ΔDβ) - Kα(ΔDα)`, `ΔFβ =
Jtot(ΔDα,ΔDβ) - Kβ(ΔDβ)` -- the α-β coupling lives entirely in `Jtot`.
"""
function _response_fock_from_densities!(ΔFα, ΔFβ, ΔDα, ΔDβ, ints)
    Jα = zeros(size(ΔFα)); Jβ = zeros(size(ΔFβ))
    Kα = zeros(size(ΔFα)); Kβ = zeros(size(ΔFβ))
    _calcJK_uhf!(Jα, Jβ, Kα, Kβ, ΔDα, ΔDβ, ints)
    Jtot = Jα .+ Jβ
    ΔFα .= Jtot .- Kα
    ΔFβ .= Jtot .- Kβ
    return ΔFα, ΔFβ
end

# Scratch-buffer core: scr.Jα/Jβ/Kα/Kβ reused instead of allocated fresh
# (and no separate Jtot temporary -- ΔFα/ΔFβ are each computed via one
# fused broadcast).
function _response_fock_from_densities!(ΔFα, ΔFβ, ΔDα, ΔDβ, ints, scr::UCPHFScratch)
    _calcJK_uhf!(scr.Jα, scr.Jβ, scr.Kα, scr.Kβ, ΔDα, ΔDβ, ints)
    ΔFα .= scr.Jα .+ scr.Jβ .- scr.Kα
    ΔFβ .= scr.Jα .+ scr.Jβ .- scr.Kβ
    return ΔFα, ΔFβ
end

"""
    build_response_fock!(ΔFα, ΔFβ, Uα, Uβ, Coα, Cvα, Coβ, Cvβ, ints)

UCPHF analog of `RHF/Hessians/CPHF.jl`'s `build_response_fock!`: builds the
symmetric trial AO densities `ΔDα = Cvα*Uα*Coα' + Coα*Uα'*Cvα'` (and the β
analog), then `_response_fock_from_densities!` for the coupled J/K
contraction.
"""
function build_response_fock!(ΔFα, ΔFβ, Uα, Uβ, Coα, Cvα, Coβ, Cvβ, ints)
    ΔDα = Cvα * Uα * Coα'
    ΔDα .+= ΔDα'
    ΔDβ = Cvβ * Uβ * Coβ'
    ΔDβ .+= ΔDβ'
    return _response_fock_from_densities!(ΔFα, ΔFβ, ΔDα, ΔDβ, ints)
end

# Scratch-buffer core: ΔDα/ΔDβ built via mul! into scr.ΔDα/scr.ΔDβ
# (through scr.tmpα/scr.tmpβ for the Cvσ*Uσ step) instead of allocating
# fresh matrix-multiply temporaries.
function build_response_fock!(ΔFα, ΔFβ, Uα, Uβ, Coα, Cvα, Coβ, Cvβ, ints, scr::UCPHFScratch)
    mul!(scr.tmpα, Cvα, Uα)
    mul!(scr.ΔDα, scr.tmpα, Coα')
    scr.ΔDα .+= scr.ΔDα'

    mul!(scr.tmpβ, Cvβ, Uβ)
    mul!(scr.ΔDβ, scr.tmpβ, Coβ')
    scr.ΔDβ .+= scr.ΔDβ'

    return _response_fock_from_densities!(ΔFα, ΔFβ, scr.ΔDα, scr.ΔDβ, ints, scr)
end

"""
    ucphf_Amatvec((Uα,Uβ), Coα, Cvα, Coβ, Cvβ, ints)

The coupled UCPHF matrix's action on a trial rotation pair, returned as a
`(nvirα,noccα)`,`(nvirβ,noccβ)` tuple matching `(Uα,Uβ)`'s shape. Thin
wrapper around `build_response_fock!` for use as a `KrylovKit.linsolve`
linear map.
"""
function ucphf_Amatvec((Uα, Uβ), Coα, Cvα, Coβ, Cvβ, ints)
    # (Uα,Uβ) as a plain Tuple works directly with KrylovKit.linsolve --
    # KrylovKit/VectorInterface handle a homogeneous Tuple of arrays as a
    # vector type natively (verified directly before relying on it here),
    # no flattening or custom vector-space glue needed.
    nbas = size(Coα, 1)
    ΔFα = zeros(nbas, nbas)
    ΔFβ = zeros(nbas, nbas)
    build_response_fock!(ΔFα, ΔFβ, Uα, Uβ, Coα, Cvα, Coβ, Cvβ, ints)
    return (Cvα' * ΔFα * Coα, Cvβ' * ΔFβ * Coβ)
end

# Scratch-buffer core of ucphf_Amatvec, called once per CG iteration --
# resultα/resultβ written in place and returned; scr.tmp2α/scr.tmp2β hold
# the Cvσ'*ΔFσ intermediate.
function ucphf_Amatvec!((resultα, resultβ), (Uα, Uβ), Coα, Cvα, Coβ, Cvβ, ints, scr::UCPHFScratch)
    build_response_fock!(scr.ΔFα, scr.ΔFβ, Uα, Uβ, Coα, Cvα, Coβ, Cvβ, ints, scr)
    mul!(scr.tmp2α, Cvα', scr.ΔFα)
    mul!(resultα, scr.tmp2α, Coα)
    mul!(scr.tmp2β, Cvβ', scr.ΔFβ)
    mul!(resultβ, scr.tmp2β, Coβ)
    return (resultα, resultβ)
end

"""
    eri_grad_JK_uhf(bset, Dtot, Dα, Dβ, iA; ij_vals=nothing, σvals=nothing)

UHF analog of `CPHF.jl`'s `eri_grad_JK`: `Jq[m,n] = Σ_rs Dtot[r,s]*d(mn|rs)/dA_q`
(Coulomb, from the total density -- shared between spins), `Kqα[m,n] =
Σ_rs Dα[r,s]*d(mr|ns)/dA_q`, `Kqβ` likewise from `Dβ`.
"""
function eri_grad_JK_uhf(bset::BasisSet, Dtot::AbstractMatrix, Dα::AbstractMatrix, Dβ::AbstractMatrix, iA::Int; ij_vals = nothing, σvals = nothing)
    # Computes the compressed derivative-quartet list
    # (GaussianBasis.∇sparseERI_2e4c) only ONCE and accumulates all three
    # (Jq,Kqα,Kqβ) in the same pass, rather than calling eri_grad_JK three
    # times (once per density) and rebuilding that list redundantly --
    # same "share the quartet list, vary only the density contraction"
    # idea RHFgrad->UHFgrad's contract_canonical/contract_cross
    # generalization uses one derivative order down.
    nbas = bset.nbas
    Jq = zeros(nbas, nbas, 3)
    Kqα = zeros(nbas, nbas, 3)
    Kqβ = zeros(nbas, nbas, 3)
    idx, xyz... = GaussianBasis.∇sparseERI_2e4c(bset, iA; ij_vals=ij_vals, σvals=σvals)
    for i in eachindex(idx)
        μ, ν, λ, σ = Int.(idx[i])
        orbit = unique([(μ, ν, λ, σ), (ν, μ, λ, σ), (μ, ν, σ, λ), (ν, μ, σ, λ),
                        (λ, σ, μ, ν), (σ, λ, μ, ν), (λ, σ, ν, μ), (σ, λ, ν, μ)])
        for q in 1:3
            ∇k = xyz[q][i]
            abs(∇k) < 1e-12 && continue
            for (a, b, c, d) in orbit
                Jq[a, b, q] += Dtot[c, d] * ∇k
                Kqα[a, c, q] += Dα[b, d] * ∇k
                Kqβ[a, c, q] += Dβ[b, d] * ∇k
            end
        end
    end
    return Jq, Kqα, Kqβ
end

"""
    ucphf_rhs(wfn::UHF, ints, iA; ij_vals=nothing, σvals=nothing)

UCPHF right-hand side, `(Bα,Bβ)` for all three Cartesian directions of atom
`iA`, alongside the `∂H,∂S,Jq_tot,Kqα,Kqβ` intermediates it was built from.
"""
function ucphf_rhs(wfn::UHF, ints, iA::Int; ij_vals = nothing, σvals = nothing)
    # ∂H/∂S/Jq_tot/Kqα/Kqβ are returned (not just B) so _uhf_hess_response
    # can reuse them instead of a second, redundant pass -- same design as
    # CPHF.jl's cphf_rhs (see ucphf_solve_full's comment). Bα/Bβ's
    # "skeleton" piece (Fskelα=∂H+Jq_tot-Kqα, Fskelβ likewise with Kqβ) and
    # occ-occ-orthonormality-correction piece (dD_S,σ:=-Coσ*S_ooσ^(y)*Coσ',
    # response Fock via _response_fock_from_densities! -- built per spin
    # but still Coulomb-coupled in its Fock response) both mirror
    # cphf_rhs's RHF formula one spin channel at a time.
    bset = BasisSet(wfn.orbitals.basis, wfn.molecule.atoms)
    nbas = bset.nbas
    Nα = wfn.molecule.Nα
    Nβ = wfn.molecule.Nβ
    nvirα = nbas - Nα
    nvirβ = nbas - Nβ
    Cα = wfn.orbitals.Cα
    Cβ = wfn.orbitals.Cβ
    Coα = @view Cα[:, 1:Nα]; Cvα = @view Cα[:, Nα+1:end]
    Coβ = @view Cβ[:, 1:Nβ]; Cvβ = @view Cβ[:, Nβ+1:end]
    epsα_o = wfn.orbitals.epsα[1:Nα]
    epsβ_o = wfn.orbitals.epsβ[1:Nβ]
    Dα = Coα * Coα'
    Dβ = Coβ * Coβ'
    Dtot = Dα .+ Dβ

    ∂H = zeros(nbas, nbas, 3)
    ∂Vtmp = zeros(nbas, nbas, 3)
    GaussianBasis.∇kinetic!(∂H, bset, iA)
    GaussianBasis.∇nuclear!(∂Vtmp, bset, iA)
    ∂H .+= ∂Vtmp

    ∂S = zeros(nbas, nbas, 3)
    GaussianBasis.∇overlap!(∂S, bset, iA)

    Jq_tot, Kqα, Kqβ = eri_grad_JK_uhf(bset, Dtot, Dα, Dβ, iA; ij_vals=ij_vals, σvals=σvals)

    Bα = zeros(nvirα, Nα, 3)
    Bβ = zeros(nvirβ, Nβ, 3)
    Gα = zeros(nbas, nbas)
    Gβ = zeros(nbas, nbas)
    for q in 1:3
        # NOTE: no factor of 2 on Jq_tot here, unlike RHF's CPHF.jl -- RHF's
        # Fskel matches F=H+2J(D)-K(D) (built from D, not the total
        # density), but UHF's Fα=H+Jtot-Kα already has Jtot built directly
        # from Dtot (eri_grad_JK_uhf's Jq IS Jq_tot already), so no extra
        # factor belongs here. (Caught via the closed-shell-vs-RHF check --
        # this is exactly the kind of bug that check exists to catch.)
        Fskelα = @view(∂H[:, :, q]) .+ (@view Jq_tot[:, :, q]) .- (@view Kqα[:, :, q])
        Fskelβ = @view(∂H[:, :, q]) .+ (@view Jq_tot[:, :, q]) .- (@view Kqβ[:, :, q])

        Sq = @view ∂S[:, :, q]
        S_oo_α = Coα' * Sq * Coα
        S_oo_β = Coβ' * Sq * Coβ
        dD_S_α = -Coα * S_oo_α * Coα'
        dD_S_β = -Coβ * S_oo_β * Coβ'
        _response_fock_from_densities!(Gα, Gβ, dD_S_α, dD_S_β, ints)

        Bα[:, :, q] .= Cvα' * Fskelα * Coα .- (Cvα' * Sq * Coα) * Diagonal(epsα_o) .+ Cvα' * Gα * Coα
        Bβ[:, :, q] .= Cvβ' * Fskelβ * Coβ .- (Cvβ' * Sq * Coβ) * Diagonal(epsβ_o) .+ Cvβ' * Gβ * Coβ
    end

    return (Bα, Bβ), ∂H, ∂S, Jq_tot, Kqα, Kqβ
end

"""
    ucphf_solve(wfn::UHF, ints, iA)

Solve the UCPHF equations for `(Uα,Uβ)`, the occupied-virtual orbital
response to displacing atom `iA`, for all three Cartesian directions. Thin
wrapper around `ucphf_solve_full` for callers that only need `(Uα,Uβ)`.
"""
function ucphf_solve(wfn::UHF, ints, iA::Int; ij_vals = nothing, σvals = nothing)
    (Uα, Uβ), _, _, _, _, _ = ucphf_solve_full(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)
    return (Uα, Uβ)
end

"""
    ucphf_solve_full(wfn::UHF, ints, iA)

Same solve as `ucphf_solve`, but also returns `ucphf_rhs`'s intermediate
`∂H,∂S,Jq_tot,Kqα,Kqβ` alongside `(Uα,Uβ)`. Solves the coupled system via
`KrylovKit.linsolve` with (per-spin Jacobi-preconditioned) CG.
"""
function ucphf_solve_full(wfn::UHF, ints, iA::Int; ij_vals = nothing, σvals = nothing)
    # The block operator [[Δεα+A^αα, A^αβ],[A^βα, Δεβ+A^ββ]] is SPD at a
    # true UHF minimum (same argument as the RHF case, one spin channel at
    # a time). Preconditioned per spin (dα:=sqrt.(εα_a-εα_i), dβ likewise)
    # from the start rather than retrofitted -- the RHF CPHF work found
    # plain CG hitting cphf_max_iter on ~85% of right-hand sides for
    # exactly this reason (the ε_a-ε_i denominators span a huge range), and
    # the same operator-conditioning issue applies per spin channel here.
    # See CPHF.jl's cphf_solve_full comment for the full derivation.
    Nα = wfn.molecule.Nα
    Nβ = wfn.molecule.Nβ
    nbas = size(wfn.orbitals.Cα, 1)
    nvirα = nbas - Nα
    nvirβ = nbas - Nβ
    Cα = wfn.orbitals.Cα
    Cβ = wfn.orbitals.Cβ
    Coα = @view Cα[:, 1:Nα]; Cvα = @view Cα[:, Nα+1:end]
    Coβ = @view Cβ[:, 1:Nβ]; Cvβ = @view Cβ[:, Nβ+1:end]
    epsα = wfn.orbitals.epsα
    epsβ = wfn.orbitals.epsβ
    Δεα = epsα[Nα+1:end] .- epsα[1:Nα]'
    Δεβ = epsβ[Nβ+1:end] .- epsβ[1:Nβ]'
    dα = sqrt.(Δεα)
    dβ = sqrt.(Δεβ)

    maxiter = Options.get("cphf_max_iter")
    tol = Options.get("cphf_conv")

    (Bα, Bβ), ∂H, ∂S, Jq_tot, Kqα, Kqβ = ucphf_rhs(wfn, ints, iA; ij_vals=ij_vals, σvals=σvals)

    Uα = zeros(nvirα, Nα, 3)
    Uβ = zeros(nvirβ, Nβ, 3)

    # Reused across every CG iteration (all 3 directions) instead of
    # letting ucphf_Amatvec/build_response_fock! allocate fresh on every
    # matvec call -- see UCPHFScratch's docstring (UHF analog of
    # CPHF.jl's cphf_Amatvec! fix).
    scr = UCPHFScratch(nbas, nvirα, Nα, nvirβ, Nβ)
    resultα = zeros(nvirα, Nα)
    resultβ = zeros(nvirβ, Nβ)

    for q in 1:3
        function precond_matvec((yα, yβ))
            ucphf_Amatvec!((resultα, resultβ), (yα ./ dα, yβ ./ dβ), Coα, Cvα, Coβ, Cvβ, ints, scr)
            return (yα .+ resultα ./ dα, yβ .+ resultβ ./ dβ)
        end
        rhs = (-Bα[:, :, q] ./ dα, -Bβ[:, :, q] ./ dβ)
        x0 = (zeros(nvirα, Nα), zeros(nvirβ, Nβ))
        sol, info = KrylovKit.linsolve(precond_matvec, rhs, x0, KrylovKit.CG(; maxiter=maxiter, tol=tol))
        Uα[:, :, q] .= sol[1] ./ dα
        Uβ[:, :, q] .= sol[2] ./ dβ
    end
    return (Uα, Uβ), ∂H, ∂S, Jq_tot, Kqα, Kqβ
end
