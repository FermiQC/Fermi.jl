"""
    Fermi.CoupledCluster.RCCD

Fermi struct that holds information about RCCSD wavefunctions

_struct tree:_

**RCCSD** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
mutable struct RCCSD{T} <: AbstractCCWavefunction
    CorrelationEnergy::T
    T1::_T where _T <: Fermi.AbstractTensor
    T2::_T where _T <: Fermi.AbstractTensor
end

# Most general function. It defines the precision and call a precision-specific function
function RCCSD()
    precision_selector = Dict{Any,Any}("single" => Float32,
                                       "double" => Float64)
    prec = precision_selector[Fermi.CurrentOptions["precision"]]
    RCCSD{prec}()
end

function RCCSD{T}() where T <: AbstractFloat
    alg_selector = Dict{Any,Any}("dpd" => DPD(),
                                 "ctf" => CTF())

    alg = alg_selector[Fermi.CurrentOptions["cc_alg"]]
    RCCD{T}(alg)
end

#implementations
include("CTF.jl")
