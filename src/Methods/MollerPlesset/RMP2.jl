struct RMP2{T} <: AbstractMPWavefunction# where T <: AbstractFloat
    CorrelationEnergy::T
    T2::_T2 where _T2 <: Fermi.AbstractTensor
end

function RMP2()
    refwfn = Fermi.HartreeFock.RHF()
    RMP2(refwfn)
end
function RMP2(ref::Fermi.HartreeFock.RHF)
    precision_selector = Dict{Any,Any}("single" => Float32,
                                       "double" => Float64)
    prec = precision_selector[Fermi.CurrentOptions["precision"]]
    RMP2{prec}(ref)
end

function RMP2{T}(ref::Fermi.HartreeFock.RHF) where T <: AbstractFloat
    # convert some options -> singletons
    alg_selector = Dict{Any,Any}("direct" => Direct(),
                                 "df" => DF(),
                                 "conv" => Conventional())
    alg = alg_selector[Fermi.CurrentOptions["mp2_type"]]
    RMP2{T}(ref,alg)
end

include("DF-RMP2.jl")
include("ConventionalRMP2.jl")
#"""
#This is the most generic constructor, so we'll put the docstring here
#"""
#function RMP2{T}(ref::Fermi.HartreeFock.RHF,alg::T2) where { T <: AbstractFloat,
#                                                                             T2 <: MP2Algorithm }
#end

# Concrete implementations of RMP2
#include("DirectRMP2.jl")
