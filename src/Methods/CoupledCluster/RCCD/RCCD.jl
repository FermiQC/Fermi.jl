mutable struct RCCD{T} <: AbstractCCWavefunction
    CorrelationEnergy::T
    MP2Energy::T
    T2::_T where _T <: Fermi.AbstractTensor
end
    

function RCCD(ref::Fermi.HartreeFock.RHF.RHFWavefunction)
    precision_selector = Dict{Any,Any}("single" => Float32,
                                       "double" => Float64)
    prec = precision_selector[Fermi.CurrentOptions["precision"]]
    RCCD{prec}(ref)
end

function RCCD{T}(ref::Fermi.HartreeFock.RHF.RHFWavefunction)
    alg_selector = Dict{Any,Any}("dpd" => DPD(),
                                 "ctf" => CTF())
    int_selector = Dict{Any,Any}("df" => DF(),
                                 "conv" => Conventional(),
                                 "cd" => CD())
    alg = alg_selector[Fermi.CurrentOptions["cc_alg"]]
    int = int_selector[Fermi.CurrentOptions["ints_type"]]
    RCCD{T}(ref,alg,int)
end

#implementations
include("RCCD_DPD_Conventional.jl")

"""
Docstring
"""
function RCCD{T}(ref::Fermi.HartreeFock.RHF.RHFWavefunction,alg::T2,int::T3) where { T <: AbstractFloat,
                                                                                    T2 <: CCAlgorithm,
                                                                                    T3 <: CCIntegralType }
end
