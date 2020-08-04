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

struct PackedTensor{T} <: AbstractTensor where T <: Number
    symmetries::Dict{Tuple,String}
    data::Array{T,1}
    ndim::Int8
    dims::Array{UInt16,1}
end

struct MemTensor{T} <: AbstractTensor where T <: Number
    data::Array{T}
    ndim::Int8
    dims::Array{UInt16,1}
end

struct DiskTensor{T} <: AbstractTensor where T <: Number
    fname::String
    buf::Array{T}
    ndim::Int8
    dims::Array{UInt32,1}
end

struct DistributedTensor{T} <: AbstractTensor where T <: Number
    data::DistributedArrays.DArray{T}
    ndim::Int8
    dims::Array{UInt32,1}
end

## MemTensor ##
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
## End MemTensor ## 
