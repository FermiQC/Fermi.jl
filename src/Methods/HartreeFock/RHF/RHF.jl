using Fermi.DIIS

using TensorOperations
using LinearAlgebra
using Formatting
import Base: show

export RHF

"""
    Fermi.HartreeFock.RHFAlgorithm

Abstract type for RHF implementations.
"""
abstract type RHFAlgorithm end


"""
    Fermi.HartreeFock.RHF

Wave function object for Restricted Hartree-Fock methods

# High Level Interface 
Run a RHF computation and return the RHF object:
```
julia> @energy rhf
```
Equivalent to
```
julia> Fermi.HartreeFock.RHF()
```
Computes RHF using information from Fermi.CurrentOptions.

# Fields

| Name   |   Description     |
|--------|---------------------|
| `molecule` |   Molecule object |
| `energy`   |   RHF Energy      |
| `ndocc`    | Number of doubly occupied spatial orbitals |
| `nvir`  | Number of virtual spatial orbitals |
| `orbitals` |    RHF Orbitals object      |
| `e_conv`   | ΔE from the last iteration  |
| `d_conv`   |  RMS from the last iteration|

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `scf_alg`      | Picks SCF algorithm               | `Int`     | [1]                   |
| `scf_max_rms`  | RMS density convergence criterion | `Float64` | [10^-9]               |
| `scf_max_iter` | Max number of iterations          | `Int`     | [50]                  |
| `scf_e_conv`   | Energy convergence criterion      | `Float64` | [10^-10]              |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `df`           | Whether to use density fitting    | `Bool`    | `true` [`false`]      |
| `jkfit`        | What aux. basis set to use for JK | `String`  | ["auto"]              |
| `diis`         | Whether to use DIIS               | `Bool`    | [`true`] `false`      |
| `oda`          | Whether to use ODA                | `Bool`    | [`true`] `false`      |
| `oda_cutoff`   | When to turn ODA off (RMS)        | `Float64` | [1E-1]                |
| `oda_shutoff`  | When to turn ODA off (iter)       | `Int`     | [20]                  |
| `scf_guess`    | Which guess density to use        | `String`  | "core" ["gwh"]        |

# Struct tree

**RHF** <: AbstractHFWavefunction <: AbstractWavefunction
"""
struct RHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    orbitals::RHFOrbitals
    e_conv::Float64
    d_conv::Float64
end

function RHF(x...)
    if !any(i-> i isa RHFAlgorithm, x)
        RHF(x..., get_scf_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for RHF method: $args"))
    end
end

# Actual RHF routine is in here
# For each implementation a singleton type must be create
struct RHFa <: RHFAlgorithm end
include("RHFa.jl")
include("RHFHelper.jl")

### MISCELLANEOUS
# Pretty printing
function string_repr(X::RHF)
    out = ""
    out = out*" ⇒ Fermi Restricted Hartree--Fock Wave function\n"
    out = out*" ⋅ Basis:                  $(X.orbitals.basis)\n"
    out = out*" ⋅ Energy:                 $(X.energy)\n"
    out = out*" ⋅ Occ. Spatial Orbitals:  $(X.ndocc)\n"
    out = out*" ⋅ Vir. Spatial Orbitals:  $(X.nvir)\n"
    out = out*"Convergence: " 
    out = out*"ΔE => $(format("{:1.2e}",abs(X.e_conv)))"
    out = out*" Dᵣₘₛ => $(format("{:1.2e}",abs(X.d_conv)))"
    return out
end

function show(io::IO, ::MIME"text/plain", X::RHF)
    print(io, string_repr(X))
end

