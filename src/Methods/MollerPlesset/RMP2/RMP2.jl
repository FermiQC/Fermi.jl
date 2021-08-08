using Fermi.HartreeFock

import Base: show

export RMP2

"""
    Fermi.HartreeFock.RHFAlgorithm

Abstract type for RMP2 implementations.
"""
abstract type RMP2Algorithm end

"""
    Fermi.MollerPlesset.get_rmp2_alg()

Returns a singleton type corresponding to a RMP2 implementation based on the options.
"""
function get_rmp2_alg()
    implemented = [RMP2a()]
    N = Options.get("mp2_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(FermiException("implementation number $N not available for RMP2."))
    end
end

"""
    RMP2{T} <: AbstractMPWavefunction

    TODO
"""
struct RMP2{T} <: AbstractMPWavefunction
    correlation::T
    energy::AbstractFloat
end

# For each implementation a singleton type must be create
struct RMP2a <: RMP2Algorithm end
include("RMP2a.jl")

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
