mutable struct RMP2{T} <: AbstractMPWavefunction where T <: AbstractFloat
    CorrelationEnergy::T
    D::T where T <: Fermi.AbstractTensor
    T2::T where T <: Fermi.AbstractTensor
end

abstract type MP2Algorithm end

struct Conventional <: MP2Algorithm end
struct DFMP2        <: MP2Algorithm end
struct Direct       <: MP2Algorithm end

function RMP2(ref::Fermi.RHF)
    RMP2{Float64}(ref)
end

function RMP2{T}(ref::Fermi.RHF) where T <: AbstractFloat
    # convert some options -> singletons
    alg_selector = Dict{Any,Any}("direct" => Direct(),
                                 "DF" => DFMP2(),
                                 "conv" => Conventional())
    alg = alg_selector[Fermi.CurrentOptions["mp2_type"]]
    RMP2{T}(ref,alg)
end

"""
This is the most generic constructor, so we'll put the docstring here
"""
function RMP2{T}(ref::Fermi.RHF,alg::T2) where { T <: AbstractFloat,
                                                T2 <: MP2Algorithm }
end

# Concrete implementations of RMP2
include("ConventionalRMP2.jl")
include("DirectRMP2.jl")
include("DF-RMP2.jl")
