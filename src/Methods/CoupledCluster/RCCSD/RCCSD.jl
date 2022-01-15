using Fermi.HartreeFock

import Base: show

export RCCSD

"""
    Fermi.CoupledCluster.RCCSDAlgorithm

Abstract type for RCCSD implementations.
"""
abstract type RCCSDAlgorithm end

"""
    Fermi.CoupledCluster.RCCSD

Wave function object for Restricted Coupled Cluster Singles and Doubles.

# High Level Interface 
Run a RCCSD computation and return the RCCSD object:
```
julia> @energy rccsd 
```
Equivalent to
```
julia> Fermi.CoupledCluster.RCCSD()
```
This function calls a constructor that runs a RCCSD computation based on the options found in `Fermi.Options`.

# Fields

| Name   |   Description     |
|--------|---------------------|
| `guessenergy` |   Energy recovered in the first iteration, normally RMP2 |
| `correlation` |   Computed RCCSD correlation energy |
| `energy`   |   Total wave function energy (Reference energy + Correlation energy)      |
| `e_conv`   | ΔE from the last iteration  |
| `t_conv`   |  Amplitudes RMS change from the last iteration|

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `cc_alg`      | Picks RCCSD algorithm              | `Int`     | [1]                   |
| `cc_e_conv`   | Energy convergence criterion           | `Float64` | [10^-10]              |
| `cc_max_rms`    | Amplitudes RMS convergence criterion   | `Float64` | [10^-10]              |
| `cc_max_iter`   | Max number of CC iterations   | `Int` | [50]              |
| `cc_damp_ratio` | Fraction of old amplitudes to be kept   | `Float64` | 0.0--1.0 [0.0]              |
| `cc_diis` | Whether to use DIIS   | `Bool` | `false` [`true`]              |
| `diis_start` | Iteration number where DIIS starts   | `Int` | [3]              |
| `cc_diis_relax` | Interval between DIIS extrapolations   | `Int` | [3]              |
| `cc_ndiis` | Maximum number of stored vectors for DIIS   | `Int` | [3]              |
| `basis`       | What basis set to use             | `String`  | ["sto-3g"]            |
| `df`          | Whether to use density fitting    | `Bool`    | `true` [`false`]      |
| `rifit`       | What aux. basis set to use for RI | `String`  | ["auto"]              |
| `drop_occ`    | Number of occupied electrons to be dropped | `Int`  | [0]              |
| `drop_vir`    | Number of virtual electrons to be dropped | `Int`  | [0]              |
"""
struct RCCSD{T} <: AbstractCCWavefunction 
    guessenergy::Float64
    correlation::T
    energy::Float64
    T1::AbstractArray{T,2}
    T2::AbstractArray{T,4}
    e_conv::T
    t_conv::T
end

"""
    Fermi.CoupledCluster.get_rccsd_alg

Returns a singleton type corresponding to a RCCSD implementation based on the options.
"""
function get_rccsd_alg(N::Int = Options.get("cc_alg"))
    try 
        return get_rccsd_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for RCCSD."))
    end
end

# For each implementation a singleton type must be create
struct RCCSDa <: RCCSDAlgorithm end
include("RCCSDa.jl")
# And a number is assigned to the implementation
get_rccsd_alg(x::Val{1}) = RCCSDa()

function RCCSD(x...)
    if !any(i-> i isa RCCSDAlgorithm, x)
        RCCSD(x..., get_rccsd_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for RCCSD method: $args"))
    end
end

## MISCELLANEOUS
# Pretty printing
function string_repr(X::RCCSD)
    out = ""
    out = out*" ⇒ Fermi Restricted CCSD Wave function\n"
    out = out*" ⋅ Correlation Energy:     $(X.correlation)\n"
    out = out*" ⋅ Total Energy:           $(X.energy)"
    return out
end

function show(io::IO, ::MIME"text/plain", X::RCCSD)
    print(io, string_repr(X))
end
