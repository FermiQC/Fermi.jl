## AbstractTensor ##
abstract type AbstractTensor end
## Methods to be implemented (see MemTensor{T} for an example
##      size
#       data
#       permutedims
#       getindex
#       zeroed initialization without intermediate array (defeats the purpose of e.g. distributed memory)

import Base.size
import Base.permutedims
import Base.getindex

## TODO: Implement zeroed initialization for all AbstractTensor

"""
    PackedTensor{T} <: AbstractTensor where T <: Number
    
Symmetry packed tesor. Not functional."""
struct PackedTensor{T} <: AbstractTensor where T <: Number
    symmetries::Dict{Tuple,String}
    data::Array{T,1}
    ndim::Int8
    dims::Array{UInt16,1}
end

## start GeneratedTensor ##
"""
    GeneratedTensor{T} <: AbstractTensor where T <: Number

Tensor for which elements are computed on the fly by `generator`. 
"""
struct GeneratedTensor{T} <: AbstractTensor where T <: Number
    generator::Function
    data::Dict{String,Any}
    ndim::Int8
    dims::Array{UInt16,1}
end

function getindex(g::GeneratedTensor,I::Vararg{Int,N}) where N
    try
        g.generator(g,I...)
    catch
        error("Incorrect number of arguments passed to getindex(::GeneratedTensor) 😦")
    end
end
## end GeneratedTensor ##

## start DiskTensor ##
struct DiskTensor{T} <: AbstractTensor where T <: Number
    fname::String
    buf::Array{T}
    ndim::Int8
    dims::Array{UInt32,1}
end
## end DiskTensor ##

struct DistributedTensor{T} <: AbstractTensor where T <: Number
    data::DistributedArrays.DArray{T}
    ndim::Int8
    dims::Array{UInt32,1}
end

## start MemTensor ##
"""
    MemTensor{T} <: AbstractTensor where T <: Number

Tensor object held entirely in memory. Thin wrapper around a standard Julia array,
but tagged as a MemTensor to fit into Fermi's dispatch system.
"""
struct MemTensor{T} <: AbstractTensor where T <: Number
    data::Array{T}
    ndim::Int8
    dims::Array{UInt16,1}
end

function MemTensor(data::Array)
    @assert eltype(data) <: Number "MemTensor was passed an array of non-numeric values."
    MemTensor{eltype(data)}(data)
end
function MemTensor{T}(data::Array{T}) where T <: Number
    ndim = convert(Int8,length(size(data)))
    dims = zeros(UInt16,ndim)
    for (idx,val) in enumerate(size(data))
        dims[idx] = convert(UInt16,val)
    end
    MemTensor{T}(data,ndim,dims)
end

function size(A::MemTensor)
    size(A.data)
end
function data(A::MemTensor)
    A.data
end

function data(A::Array)
    A
end

function permutedims(A::MemTensor,tup)
    A.data .= permutedims(A.data,tup)
end

function getindex(A::MemTensor,I...)
    return A.data[I...]
end
## end MemTensor ## 
