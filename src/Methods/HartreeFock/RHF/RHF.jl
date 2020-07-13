abstract type RHFAlgorithm end
struct ConventionalRHF <: RHFAlgorithm end

mutable struct RHFWavefunction{T} <: AbstractHFWavefunction where T <: AbstractFloat
    molecule::Molecule
    energy::T
    nocc::Int
    nvir::Int
    C::Array{T,2} 
    eps::Array{T,1} 
end

include("ConventionalRHF.jl")
