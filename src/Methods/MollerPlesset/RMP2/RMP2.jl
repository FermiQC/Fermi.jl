using Fermi.HartreeFock

import Base: show

export RMP2

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

include("RMP2a.jl")

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
