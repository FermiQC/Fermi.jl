abstract type RHFAlgorithm end
struct ConventionalRHF <: RHFAlgorithm end

mutable struct RHF{T} <: AbstractHFWavefunction where T <: AbstractFloat
    molecule::Molecule
    energy::T
    nocc::Int
    nvir::Int
    C::Array{T,2} 
    F::Array{T,2}
end

include("ConventionalRHF.jl")
