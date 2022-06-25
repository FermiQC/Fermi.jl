using Fermi.HartreeFock

import Base: show

export RMP2

"""
    Fermi.MollerPlesset.RMP2Algorithm

Abstract type for RMP2 implementations.
"""
abstract type RMP2Algorithm end

"""
    Fermi.MollerPlesset.RMP2

Wave function object for Second Order Restricted Moller-Plesset methods.

# High Level Interface 
Run a RMP2 computation and return the RMP2 object:
```
julia> @energy rmp2
```
Equivalent to
```
julia> Fermi.MollerPlesset.RMP2()
```
This function calls a constructor that runs a MP2 computation based on the options found in `Fermi.Options`.

# Fields

| Name   |   Description     |
|--------|---------------------|
| `correlation` |   Computed RMP2 correlation energy |
| `energy`   |   Total wave function energy (Reference energy + Correlation energy)      |

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `mp2_alg`      | Picks MP2 algorithm               | `Int`     | [1]                   |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `df`           | Whether to use density fitting    | `Bool`    | `true` [`false`]      |
| `rifit`        | What aux. basis set to use for RI | `String`  | ["auto"]              |
| `drop_occ`    | Number of occupied electrons to be dropped | `Int`  | [0]              |
| `drop_vir`    | Number of virtual electrons to be dropped | `Int`  | [0]              |
"""
struct RMP2{T} <: AbstractMPWavefunction
    correlation::T
    energy::T
end

"""
    Fermi.MollerPlesset.get_rmp2_alg()

Returns a singleton type corresponding to a RMP2 implementation based on the options.
"""
function get_rmp2_alg(N::Int = Options.get("mp2_alg"))
    try 
        return get_rmp2_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for RMP2."))
    end
end

# For each implementation a singleton type must be create
struct RMP2a <: RMP2Algorithm end
include("RMP2a.jl")
# And a number is assigned to the implementation
get_rmp2_alg(x::Val{1}) = RMP2a()

struct RMP2b <: RMP2Algorithm end
include("RMP2b.jl")
# And a number is assigned to the implementation
get_rmp2_alg(x::Val{2}) = RMP2b()

function RMP2(x...)
    if !any(i-> i isa RMP2Algorithm, x)
        RMP2(x..., get_rmp2_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for RMP2 method: $args"))
    end
end

# Gradient methods
include("Gradients/RMP2grad.jl")

## MISCELLANEOUS
# Pretty printing
function string_repr(X::RMP2)
    out = ""
    out = out*" ⇒ Fermi Restricted MP2 Wave function\n"
    out = out*" ⋅ Correlation Energy:     $(X.correlation)\n"
    out = out*" ⋅ Total Energy:           $(X.energy)"
    return out
end

function show(io::IO, ::MIME"text/plain", X::RMP2)
    print(io, string_repr(X))
end
