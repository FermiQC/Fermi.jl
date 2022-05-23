using Fermi

import Base: show

export FCI

"""
    Fermi.ConfigurationInteraction.FCIAlgorithm

Abstract type for FCI implementations.
"""
abstract type FCIAlgorithm end

"""
    Fermi.ConfigurationInteraction.RFCI

Wave function object for Restricted Full Configuration Interaction.

# High Level Interface 
Run a RFCI computation and return the RFCI object:
```
julia> @energy rfci
```
Equivalent to
```
julia> Fermi.ConfigurationInteraction.RFCI()
```
This function calls a constructor that runs a RFCI computation based on the options found in `Fermi.Options`.

# Fields

| Name   |   Description     |
|--------|---------------------|
| `energy`   |   Total wave function energy (Reference energy + Correlation energy)      |

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `drop_occ`    | Number of occupied electrons to be dropped | `Int`  | [0]              |
| `drop_vir`    | Number of virtual electrons to be dropped | `Int`  | [0]              |
"""
struct RFCI{T} <: AbstractCCWavefunction 
    energy::Float64
    correlation::T
end

"""
    Fermi.CoupledCluster.get_rfci_alg

Returns a singleton type corresponding to a RFCI implementation based on the options.
"""
function get_rfci_alg(N::Int = Options.get("ci_alg"))
    try 
        return get_rfci_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for RFCI."))
    end
end

# For each implementation a singleton type must be create
struct RFCIa <: FCIAlgorithm end
include("RFCI.jl")
# And a number is assigned to the implementation
get_rfci_alg(x::Val{1}) = RFCIa()

function RFCI(x...)
    if !any(i-> i isa FCIAlgorithm, x)
        RFCI(x..., get_rfci_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for RFCI method: $args"))
    end
end

## MISCELLANEOUS
# Pretty printing
function string_repr(X::RFCI)
    out = ""
    out = out*" ⇒ Fermi Restricted FCI Wave function\n"
    out = out*" ⋅ Correlation Energy:     $(X.correlation)\n"
    out = out*" ⋅ Total Energy:           $(X.energy)"
    return out
end

function show(io::IO, ::MIME"text/plain", X::RFCI)
    print(io, string_repr(X))
end