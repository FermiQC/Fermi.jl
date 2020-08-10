struct RMP3{T} <: AbstractMPWavefunction# where T <: AbstractFloat
    CorrelationEnergy::T
    D::_T2 where _T2 <: Fermi.AbstractTensor
    T2::_T2 where _T2 <: Fermi.AbstractTensor
end

include("DF-RMP3.jl")

function RMP3()
    refwfn = Fermi.HartreeFock.RHF()
    RMP3(refwfn)
end
function RMP3(ref::Fermi.HartreeFock.RHF)
    precision_selector = Dict{Any,Any}("single" => Float32,
                                       "double" => Float64)
    prec = precision_selector[Fermi.CurrentOptions["precision"]]
    RMP3{prec}(ref)
end

function RMP3{T}(ref::Fermi.HartreeFock.RHF) where T <: AbstractFloat
    # convert some options -> singletons
    alg_selector = Dict{Any,Any}("direct" => Direct(),
                                 "df" => DF(),
                                 "conv" => Conventional())
    alg = alg_selector[Fermi.CurrentOptions["mp2_type"]]
    RMP3{T}(ref,alg)
end

