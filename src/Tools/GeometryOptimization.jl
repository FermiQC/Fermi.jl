using LinearAlgebra
using Molecules
using Optim

export @optimize

# ---------------------------------------------------------------------------
# Gauge-fixing: remove the 3 translational + 3 rotational degrees of freedom
# from a Cartesian geometry by pinning three anchor atoms (a1 at the origin,
# a2 on the +x axis, a3 in the xy-plane), instead of using full internal
# coordinates. See the geometry optimizer plan/discussion for why this is a
# valid *exact* global reparametrization (not a small-displacement/tangent-space
# approximation): any non-degenerate 3-point configuration can always be
# rotated into this canonical frame, and we never need to re-derive the gauge
# mid-optimization since we optimize directly in it.
# ---------------------------------------------------------------------------

"""
    select_anchors(atoms)

Pick three well-conditioned anchor atoms (a1, a2, a3) used to gauge-fix
translation/rotation: a1 is atom 1, a2 is the atom farthest from a1 (a
well-separated axis), a3 is the atom that maximizes perpendicular distance
from the a1-a2 axis (best-conditioned for defining a plane). Returns
`(a1, a2, nothing)` for a 2-atom system (no plane needed). Throws a
`FermiException` if the molecule is linear (no valid a3 exists) -- linear
molecules have 3N-5 internal degrees of freedom, not 3N-6, and are not yet
supported by this gauge-fixing scheme.
"""
function select_anchors(atoms)
    N = length(atoms)
    N < 2 && throw(FermiException("Geometry optimization requires at least 2 atoms."))

    a1 = 1
    dists = [norm(atoms[i].xyz - atoms[a1].xyz) for i in eachindex(atoms)]
    dists[a1] = -Inf
    a2 = argmax(dists)

    N == 2 && return (a1, a2, nothing)

    e1 = normalize(atoms[a2].xyz - atoms[a1].xyz)
    perp_dists = Float64[]
    for i in eachindex(atoms)
        if i == a1 || i == a2
            push!(perp_dists, -Inf)
            continue
        end
        v = atoms[i].xyz - atoms[a1].xyz
        push!(perp_dists, norm(v - (v⋅e1)*e1))
    end
    a3 = argmax(perp_dists)

    if perp_dists[a3] < 1e-3
        throw(FermiException(
            "Molecule appears to be linear. Cartesian gauge-fixed geometry " *
            "optimization does not yet support linear molecules (they have " *
            "3N-5 internal degrees of freedom, not 3N-6). This is a known " *
            "limitation -- internal coordinates would be needed to handle it."
        ))
    end

    return (a1, a2, a3)
end

"""
    align_to_gauge(atoms, a1, a2, a3)

One-time rigid rotation+translation of `atoms` so that `a1` sits at the
origin, `a2` sits on the +x axis, and `a3` sits in the xy-plane (z=0). Returns
a new `Vector{Atom}` (same Z/mass, rotated xyz, still in Angstrom).
"""
function align_to_gauge(atoms, a1, a2, a3)
    t = atoms[a1].xyz
    e1 = normalize(atoms[a2].xyz - t)

    if a3 === nothing
        # Only 2 atoms: any e2,e3 completing the frame works, nothing else to align.
        e2 = abs(e1[1]) < 0.9 ? normalize([1.0,0,0] - (e1[1])*e1) : normalize([0,1.0,0] - (e1[2])*e1)
    else
        v3 = atoms[a3].xyz - t
        e2 = normalize(v3 - (v3⋅e1)*e1)
    end
    e3 = cross(e1, e2)
    R = permutedims(hcat(e1, e2, e3))

    return [Atom(atom.Z, atom.mass, Vector(R*(atom.xyz - t))) for atom in atoms]
end

"""
    pack(atoms_aligned, a1, a2, a3, others) -> Vector{Float64}

Extract the free coordinates (Bohr) from a geometry already in the canonical
gauge: `a2`'s x, `a3`'s x and y (if `a3 !== nothing`), then every other atom's
full x,y,z, in the order given by `others`.
"""
function pack(atoms_aligned, a1, a2, a3, others)
    x = Float64[atoms_aligned[a2].xyz[1]]
    if a3 !== nothing
        push!(x, atoms_aligned[a3].xyz[1], atoms_aligned[a3].xyz[2])
    end
    for i in others
        append!(x, atoms_aligned[i].xyz)
    end
    return x ./ Molecules.bohr_to_angstrom
end

"""
    unpack(x, atoms_meta, a1, a2, a3, others, charge, multiplicity) -> Molecule

Inverse of `pack`: rebuild a full `Molecule` (Angstrom, matching `Atom.xyz`'s
convention) from the free coordinates `x` (Bohr). `atoms_meta` supplies Z/mass
for every atom; only positions come from `x`.
"""
function unpack(x, atoms_meta, a1, a2, a3, others, charge, multiplicity)
    x_ang = x .* Molecules.bohr_to_angstrom
    N = length(atoms_meta)
    newxyz = Vector{Vector{Float64}}(undef, N)
    newxyz[a1] = [0.0, 0.0, 0.0]

    idx = 1
    newxyz[a2] = [x_ang[idx], 0.0, 0.0]; idx += 1
    if a3 !== nothing
        newxyz[a3] = [x_ang[idx], x_ang[idx+1], 0.0]; idx += 2
    end
    for i in others
        newxyz[i] = x_ang[idx:idx+2]; idx += 3
    end

    newatoms = [Atom(atoms_meta[i].Z, atoms_meta[i].mass, newxyz[i]) for i in 1:N]
    return Molecule(newatoms, charge, multiplicity)
end

"""
    pack_gradient(∂E, a1, a2, a3, others) -> Vector{Float64}

Project the full (Natoms x 3, Hartree/Bohr) Cartesian gradient onto the free
coordinates, in the same order as `pack`. `unpack` is an affine embedding
(free values inserted at fixed positions, the rest held at constants), so this
is exactly the corresponding sub-vector of the full gradient -- no extra
factors needed.
"""
function pack_gradient(∂E, a1, a2, a3, others)
    g = Float64[∂E[a2,1]]
    if a3 !== nothing
        push!(g, ∂E[a3,1], ∂E[a3,2])
    end
    for i in others
        append!(g, ∂E[i,:])
    end
    return g
end

# ---------------------------------------------------------------------------
# Method table: energy/gradient/warm-start functions per supported method.
# Only RHF and UHF for now -- warm-starting relies on the wavefunction-guess
# projection dispatches added to RHFa.jl/UHFa.jl specifically for this, which
# post-HF methods don't have yet.
# ---------------------------------------------------------------------------

function _rhf_energy(mol::Molecule)
    Fermi.HartreeFock.RHF(mol)
end
function _rhf_energy_warm(prev_wfn, mol::Molecule)
    target_ints = Fermi.Integrals.IntegralHelper{Float64}(molecule=mol)
    Fermi.HartreeFock.RHF(prev_wfn, target_ints, Fermi.HartreeFock.get_rhf_alg())
end

function _uhf_energy(mol::Molecule)
    Fermi.HartreeFock.UHF(mol)
end
function _uhf_energy_warm(prev_wfn, mol::Molecule)
    target_ints = Fermi.Integrals.IntegralHelper{Float64}(molecule=mol)
    Fermi.HartreeFock.UHF(prev_wfn, target_ints, Fermi.HartreeFock.get_uhf_alg())
end

# grad(wfn, mol): RHF has an analytic gradient that reuses the converged
# wavefunction directly (no extra SCF). UHF only has a findif gradient, which
# always restarts SCF at displaced geometries built from `mol` and has no way
# to reuse `wfn` at all -- so it ignores `wfn` and takes the extra cost that's
# inherent to finite differences.
_rhf_grad(wfn, mol) = Fermi.HartreeFock.RHFgrad(wfn)
_uhf_grad(wfn, mol) = Fermi.HartreeFock.UHFgrad(mol)

const opt_dict = Dict{String, NamedTuple}(
    "rhf" => (energy = _rhf_energy, energy_warm = _rhf_energy_warm, grad = _rhf_grad),
    "uhf" => (energy = _uhf_energy, energy_warm = _uhf_energy_warm, grad = _uhf_grad),
)

# ---------------------------------------------------------------------------
# Core driver
# ---------------------------------------------------------------------------

"""
    Fermi.optimize_geometry(method::String, mol::Molecule=Molecule())

Optimize `mol`'s geometry at the given `method` ("rhf" or "uhf"), using
Cartesian coordinates with translation/rotation gauge-fixed away (see
`select_anchors`/`align_to_gauge`) and `Optim.jl`'s BFGS for the underlying
quasi-Newton step. Returns the final converged wavefunction (at the optimized
geometry).
"""
function optimize_geometry(method::String, mol::Molecule=Molecule())
    haskey(opt_dict, method) || throw(FermiException("Invalid or unsupported method for geometry optimization: \"$method\". Implemented methods: $(join(sort(collect(keys(opt_dict))), ", "))"))
    funcs = opt_dict[method]

    a1, a2, a3 = select_anchors(mol.atoms)
    others = [i for i in eachindex(mol.atoms) if i ∉ (a1, a2, a3)]

    atoms0 = align_to_gauge(mol.atoms, a1, a2, a3)
    x0 = pack(atoms0, a1, a2, a3, others)

    maxit      = Options.get("geom_max_iter")
    e_conv     = Options.get("geom_e_conv")
    grms_conv  = Options.get("geom_grms_conv")
    gmax_conv  = Options.get("geom_gmax_conv")
    drms_conv  = Options.get("geom_drms_conv")
    dmax_conv  = Options.get("geom_dmax_conv")

    # Memoized energy+gradient evaluation: only recompute (and only rerun SCF)
    # when x actually changes, regardless of the order/pattern in which Optim.jl
    # calls the separate f(x) / g!(G,x) functions below (its line search may
    # evaluate f at several trial points per iteration; g! is called once at
    # whichever point is finally accepted). Also chains warm-starts between
    # *every* distinct geometry evaluated, not just accepted iterates.
    cache_x = Ref{Vector{Float64}}(fill(NaN, length(x0)))
    cache_wfn = Ref{Any}(nothing)

    function ensure!(x)
        if cache_wfn[] === nothing || cache_x[] != x
            newmol = unpack(x, mol.atoms, a1, a2, a3, others, mol.charge, mol.multiplicity)
            wfn = cache_wfn[] === nothing ? funcs.energy(newmol) : funcs.energy_warm(cache_wfn[], newmol)
            cache_wfn[] = wfn
            cache_x[] = copy(x)
        end
        return cache_wfn[]
    end

    fE(x) = ensure!(x).energy
    function gE!(G, x)
        wfn = ensure!(x)
        newmol = unpack(x, mol.atoms, a1, a2, a3, others, mol.charge, mol.multiplicity)
        ∂E = funcs.grad(wfn, newmol)
        G .= pack_gradient(∂E, a1, a2, a3, others)
    end

    ite = Ref(0)
    function cb(state)
        ite[] += 1
        if isnan(state.f_x_previous)
            output(" {:>5} {:>18.10f} {:>12} {:>12} {:>12}", ite[], state.f_x, "--", "--", "--")
            return false
        end

        ΔE = abs(state.f_x - state.f_x_previous)
        g = state.g_x
        grms = √(sum(abs2, g)/length(g))
        gmax = maximum(abs.(g))
        dx = state.x - state.x_previous
        drms = √(sum(abs2, dx)/length(dx))
        dmax = maximum(abs.(dx))

        output(" {:>5} {:>18.10f} {:>12.4e} {:>12.4e} {:>12.4e}", ite[], state.f_x, ΔE, grms, dmax)

        return ΔE < e_conv && grms < grms_conv && gmax < gmax_conv && drms < drms_conv && dmax < dmax_conv
    end

    output("\n Geometry optimization ({})", uppercase(method))
    output(" {:>5} {:>18} {:>12} {:>12} {:>12}", "Iter", "Energy", "ΔE", "GRMS", "Dmax")
    output(repeat("-", 65))

    # Seed the inverse Hessian with a uniform diagonal force-constant guess
    # (~0.5 Hartree/Bohr^2, a standard simple heuristic) instead of identity,
    # so the first step isn't wild before BFGS has learned any real curvature.
    result = Optim.optimize(fE, gE!, x0,
        BFGS(initial_invH = x -> Matrix(2.0*I, length(x), length(x))),
        Optim.Options(iterations = maxit, g_tol = 0.0, callback = cb))

    output(repeat("-", 65))

    final_x = Optim.minimizer(result)
    final_wfn = ensure!(final_x)
    return final_wfn
end

"""
    Fermi.@optimize

Macro to run a geometry optimization given current options, mirroring
`@energy`/`@gradient`. Convergence and iteration limits are controlled with
`@set geom_max_iter <N>`, `@set geom_e_conv <val>`, etc.

# Examples
```
wfn = @optimize rhf
```
Starting from an explicit `Molecule`:
```
wfn = @optimize mol => rhf
```

# Implemented methods: rhf, uhf
"""
macro optimize(comm)
    clean_up(s) = String(filter(c->!occursin(c," :"),s))
    A = clean_up(repr(comm))
    while A[1] == '(' && A[end] == ')'
        A = A[2:end-1]
    end

    if occursin("=>", A)
        arg, method = split(A, "=>")
    elseif occursin("<=", A)
        method, arg = split(A, "<=")
    else
        method = A
        arg = ""
    end

    method = lowercase(String(method))
    arg = String(arg)

    call = isempty(arg) ? "Fermi.optimize_geometry(\"$method\")" : "Fermi.optimize_geometry(\"$method\", $arg)"
    expr_out = Meta.parse(call)

    return esc(expr_out)
end
