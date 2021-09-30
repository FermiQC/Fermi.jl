abstract type RpTAlgorithm end

function get_rpt_alg()
    implemented = [RpTa(), RpTb()]
    N = Options.get("pt_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(FermiException("implementation number $N not available for RCCSD(T)."))
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
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for RCCSD(T) method: $args"))
    end
end

#implementations
struct RpTa <: RpTAlgorithm end
struct RpTb <: RpTAlgorithm end
include("RpTa.jl")
include("RpTb.jl")