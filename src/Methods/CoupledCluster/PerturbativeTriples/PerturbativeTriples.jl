abstract type RpTAlgorithm end

function get_rpt_alg()
    implemented = [RpTa()]
    N = Options.get("pt_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(InvalidFermiOption("implementation number $N not available for RCCSD(T)."))
    end
end
"""
    Fermi.CoupledCluster.RCCSDpT

Fermi struct that holds information about RCCSD(T) wavefunctions

_struct tree:_

**RCCSD(T)** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct RCCSDpT{T} <: AbstractCCWavefunction
    CCSD::RCCSD{T}
    energy::T
    correction::T
end

function RCCSDpT(x...)
    if !any(i-> i isa RpTAlgorithm, x)
        RCCSDpT(x..., get_rpt_alg())
    else
        throw(MethodArgument("invalid arguments for RCCSD(T) method: $(x[1:end-1])"))
    end
end

#implementations
struct RpTa <: RpTAlgorithm end
include("RpTa.jl")