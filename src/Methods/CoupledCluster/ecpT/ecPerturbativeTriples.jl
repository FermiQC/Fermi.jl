using Fermi.ConfigurationInteraction: Determinant

"""
    Fermi.CoupledCluster.RCCSDpT

Fermi struct that holds information about RCCSD(T) wavefunctions

_struct tree:_

**RCCSD(T)** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct ecRCCSDpT{T} <: AbstractCCWavefunction
    ecCCSD::ecRCCSD{T}
    correction::T
end

function select_precision(A::String)
    implemented = Dict{Any,Any}("single" => Float32,
                                "double" => Float64)
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid precision: $(A)"))
    end
end

#implementations
include("ijk.jl")
