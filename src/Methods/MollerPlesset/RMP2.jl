mutable struct RMP2{T} <: AbstractMPWavefunction where T <: AbstractFloat
    CorrelationEnergy::T
    D::_T2 where _T2 <: Fermi.AbstractTensor
    T2::_T2 where _T2 <: Fermi.AbstractTensor
end


abstract type MP2Algorithm end

struct Conventional <: MP2Algorithm end
struct DF           <: MP2Algorithm end
struct Direct       <: MP2Algorithm end


function RMP2(ref::Fermi.HartreeFock.RHF.RHFWavefunction)
    RMP2{Float64}(ref)
end

function RMP2{T}(ref::Fermi.HartreeFock.RHF.RHFWavefunction) where T <: AbstractFloat
    # convert some options -> singletons
    alg_selector = Dict{Any,Any}("direct" => Direct(),
                                 "DF" => DF(),
                                 "conv" => Conventional())
    alg = alg_selector[Fermi.CurrentOptions["mp2_type"]]
    RMP2{T}(ref,alg)
end

include("ConventionalRMP2.jl")

"""
This is the most generic constructor, so we'll put the docstring here
"""
function RMP2{T}(ref::Fermi.HartreeFock.RHF.RHFWavefunction,alg::T2) where { T <: AbstractFloat,
                                                                             T2 <: MP2Algorithm }
end

# Concrete implementations of RMP2
#include("DirectRMP2.jl")
#include("DF-RMP2.jl")
