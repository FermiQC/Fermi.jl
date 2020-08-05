using Fermi.ConfigurationInteraction: Determinant, excitation_level, αexcitation_level, βexcitation_level, αexclusive, βexclusive, phase, create, annihilate

"""
    Fermi.CoupledCluster.ecRCCSD

Fermi struct that holds information about ecRCCSD wavefunctions

# Fields

    CorrelationEnergy   CCSD correlation energy
    T1                  T1 amplitudes
    T2                  T2 amplitudes

_struct tree:_

**RCCSD** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct ecRCCSD{T} <: AbstractCCWavefunction
    GuessEnergy::T
    CorrelationEnergy::T
    T1::_T where _T <: Fermi.AbstractTensor
    T2::_T where _T <: Fermi.AbstractTensor
end

#function select_precision(A::String)
#    implemented = Dict{Any,Any}("single" => Float32,
#                                "double" => Float64)
#    try
#        return implemented[A]
#    catch KeyError
#        throw(Fermi.InvalidFermiOptions("Invalid precision: $(A)"))
#    end
#end
#
#function select_algorithm(A::String)
#    implemented = Dict{Any,Any}("DPD" => DPD(),
#                                "CTF" => CTF())
#    try
#        return implemented[A]
#    catch KeyError
#        throw(Fermi.InvalidFermiOptions("Invalid CC algorithm: $(A)"))
#    end
#end

#implementations
include("CTF.jl")

# Most general function. It defines the precision and call a precision-specific function
"""
    Fermi.CoupledCluster.RCCSD()

Compute a RCCSD wave function using Fermi.CurrentOptions data.
"""
function ecRCCSD()
    @output "Selecting precision...\n"
    prec = select_precision(Fermi.CurrentOptions["precision"])
    ecRCCSD{prec}()
end

"""
    Fermi.CoupledCluster.RCCSD{T}()

Compute a RCCSD wave function for a given precision T (Float64 or Float32)
"""
function ecRCCSD{T}() where T <: AbstractFloat
    @output "Selecting Algorithm...\n"
    alg = select_algorithm(Fermi.CurrentOptions["cc_alg"])
    ecRCCSD{T}(alg)
end
