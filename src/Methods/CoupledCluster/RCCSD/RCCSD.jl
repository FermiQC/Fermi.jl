"""
    Fermi.CoupledCluster.RCCD

Fermi struct that holds information about RCCSD wavefunctions

_struct tree:_

**RCCSD** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct RCCSD{T} <: AbstractCCWavefunction
    CorrelationEnergy::T
    T1::_T where _T <: Fermi.AbstractTensor
    T2::_T where _T <: Fermi.AbstractTensor
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

function select_algorithm(A::String)
    implemented = Dict{Any,Any}("DPD" => DPD(),
                                "CTF" => CTF())
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid CC algorithm: $(A)"))
    end
end

#implementations
include("CTF.jl")

# Most general function. It defines the precision and call a precision-specific function
function RCCSD()
    println("Selecting precision...")
    prec = select_precision(Fermi.CurrentOptions["precision"])
    RCCSD{prec}()
end

function RCCSD{T}() where T <: AbstractFloat
    println("Selecting Algorithm...")
    alg = select_algorithm(Fermi.CurrentOptions["cc_alg"])
    RCCSD{T}(alg)
end

