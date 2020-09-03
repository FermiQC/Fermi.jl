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

# Most general function. It defines the precision and call a precision-specific function
#"""
#    Fermi.CoupledCluster.RCCSD()
#
#Compute a RCCSD wave function using Fermi.CurrentOptions data.
#"""
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

# kernel (loop logic) from main RCCSD
include("../RCCSD/kernel.jl")

#implementations
include("CTF.jl")
include("ClusterDecomposition.jl")
