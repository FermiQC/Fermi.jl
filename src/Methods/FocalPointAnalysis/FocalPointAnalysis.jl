module FocalPointAnalysis

using Fermi
using Fermi.Output
using Fermi.Integrals: ConventionalAOIntegrals, PhysRestrictedMOIntegrals
using Fermi.HartreeFock: RHF
using Fermi.Geometry: Molecule
using Fermi.CoupledCluster: RCCSD, RCCSDpT
using PrettyTables

abstract type AbstractFPA end

struct FPA{T} <: AbstractFPA
    basis::Array{String,1}
    methods::Array{String,1}
    data::Array{T,2}
    extrapolated_indexes::Array{Tuple{Int64,Int64},1}
end

struct ecFPA{T} <: AbstractFPA
    basis::Array{String,1}
    methods::Array{String,1}
    data::Array{T,2}
    extrapolated_indexes::Array{Tuple{Int64,Int64},1}
end

function print_fpa(fpa::T) where T <: AbstractFPA
    format = ft_printf("%15.10f")
    put_braket = (v,i,j) -> (i,j) in fpa.extrapolated_indexes ? "["*v*"]" : v

    out = pretty_table(String, fpa.data, fpa.methods, alignment=:c,
    formatters = (format, put_braket), row_names=fpa.basis, 
    row_name_column_title="Basis Set")

    @output "\n{}\n" out
end

include("Extrapolation.jl")
include("FPA.jl")
include("ecFPA.jl")

end #module
