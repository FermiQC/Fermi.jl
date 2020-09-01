using LinearAlgebra
"""
    Fermi.CoupledCluster.RCCSD

Fermi struct that holds information about RCCSD wavefunctions

# Fields

    CorrelationEnergy   CCSD correlation energy
    T1                  T1 amplitudes
    T2                  T2 amplitudes

_struct tree:_

**RCCSD** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct BCCD{T} <: AbstractCCWavefunction
    GuessEnergy::T
    CorrelationEnergy::T
    T1::_T where _T <: Fermi.AbstractTensor
    T2::_T where _T <: Fermi.AbstractTensor
end

# Most general function. It defines the precision and call a precision-specific function
"""
    Fermi.CoupledCluster.RCCSD()

Compute a RCCSD wave function using Fermi.CurrentOptions data.
"""
function BCCD()
    prec = select_precision(Fermi.CurrentOptions["precision"])
    dummy = BCCD{prec}(0.0,0.0,Fermi.MemTensor(zeros(prec,0,0)),Fermi.MemTensor(zeros(prec,0,0,0,0)))
    BCCD(dummy)
end

function BCCD(guess::BCCD{T}) where T <: AbstractFloat
    prec = eltype(guess.T2.data)
    BCCD{prec}(guess)
end

"""
    Fermi.CoupledCluster.RCCSD{T}()

Compute a RCCSD wave function for a given precision T (Float64 or Float32)
"""
function BCCD{T}(guess::BCCD{Tb}) where { T <: AbstractFloat,
                                           Tb <: AbstractFloat }
    alg = select_algorithm(Fermi.CurrentOptions["cc_alg"])
    BCCD{T}(guess,alg)
end



# implementation specific functions (look at one of these for guidance implementing a new method!)
include("CTF.jl")
include("DF-CTF.jl")

# main kernel (ccsd iteration logic)
include("kernel.jl")
