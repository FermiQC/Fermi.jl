using Fermi.HartreeFock

import Base: show

export ecRCCSD

"""
    Fermi.CoupledCluster.ecRCCSDAlgorithm

Abstract type for ecRCCSD implementations.
"""
abstract type ecRCCSDAlgorithm end

"""
    Fermi.CoupledCluster.ecRCCSD

Wave function object for Restricted Coupled Cluster Singles and Doubles.

# High Level Interface 
Run a ecRCCSD computation and return the ecRCCSD object:
```
julia> @energy rccsd 
```
Equivalent to
```
julia> Fermi.CoupledCluster.ecRCCSD()
```
This function calls a constructor that runs a ecRCCSD computation based on the options found in `Fermi.Options`.

# Fields

| Name   |   Description     |
|--------|---------------------|
| `guessenergy` |   Energy recovered in the first iteration, normally RMP2 |
| `correlation` |   Computed ecRCCSD correlation energy |
| `energy`   |   Total wave function energy (Reference energy + Correlation energy)      |
| `e_conv`   | ΔE from the last iteration  |
| `t_conv`   |  Amplitudes RMS change from the last iteration|

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `cc_alg`      | Picks ecRCCSD algorithm              | `Int`     | [1]                   |
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
struct ecRCCSD{T} <: AbstractCCWavefunction 
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

Returns a singleton type corresponding to a ecRCCSD implementation based on the options.
"""
function get_recccsd_alg(N::Int = Options.get("cc_alg"))
    try 
        return get_recccsd_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for ecRCCSD."))
    end
end

# For each implementation a singleton type must be create
struct ecRCCSDa <: ecRCCSDAlgorithm end
include("ClusterDecomposition.jl")
include("ecRCCSDa.jl")
# And a number is assigned to the implementation
get_recccsd_alg(x::Val{1}) = ecRCCSDa()

function ecRCCSD(x...)
    if !any(i-> i isa ecRCCSDAlgorithm, x)
        ecRCCSD(x..., get_recccsd_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for ecRCCSD method: $args"))
    end
end

## MISCELLANEOUS
# Pretty printing
function string_repr(X::ecRCCSD)
    out = ""
    out = out*" ⇒ Fermi Restricted CCSD Wave function\n"
    out = out*" ⋅ Correlation Energy:     $(X.correlation)\n"
    out = out*" ⋅ Total Energy:           $(X.energy)"
    return out
end

function show(io::IO, ::MIME"text/plain", X::ecRCCSD)
    print(io, string_repr(X))
end