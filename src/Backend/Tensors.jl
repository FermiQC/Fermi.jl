abstract type AbstractTensor end
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

