module FocalPointAnalysis

using Fermi
using Fermi.Output
using Fermi.Integrals: IntegralHelper
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

function print_fpa(fpa::T; incremental=false) where T <: AbstractFPA

    data = copy(fpa.data)
    if incremental
        for k in size(fpa.data,2):-1:2
           data[:,k] -= data[:,k-1]
        end
    end

    format = ft_printf("%15.10f")
    put_braket = (v,i,j) -> (i,j) in fpa.extrapolated_indexes ? "["*v*"]" : v

    out = pretty_table(String, data, fpa.methods, alignment=:c,
    formatters = (format, put_braket), row_names=fpa.basis, 
    row_name_column_title="Basis Set")

    @output "\n{}\n" out
end

function Base.:*(N::T, F::FPA) where T <: Number
    return FPA(F.basis, F.methods, N.*F.data, F.extrapolated_indexes)
end

function Base.:+(A::FPA, B::FPA)
    @assert A.methods == B.methods
    @assert A.basis == B.basis
    @assert A.extrapolated_indexes == B.extrapolated_indexes
    return FPA(A.basis, B.methods, A.data+B.data, A.extrapolated_indexes)
end

function Base.:-(A::FPA, B::FPA)
    @assert A.methods == B.methods
    @assert A.basis == B.basis
    @assert A.extrapolated_indexes == B.extrapolated_indexes
    return FPA(A.basis, B.methods, A.data-B.data, A.extrapolated_indexes)
end

include("Extrapolation.jl")
include("FPA.jl")
include("ecFPA.jl")

end #module
