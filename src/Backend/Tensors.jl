abstract type AbstractTensor{T} end
struct MemTensor{T} <: AbstractTensor{T} where T <: Number
    data::Array{T}
    ndim::Int8
    dims::Array{UInt16,1}
end

struct DiskTensor{T} <: AbstractTensor{T} where T <: Number
    fname::String
    buf::Array{T}
    ndim::Int8
    dims::Array{UInt32,1}
end

struct DistributedTensor{T} <: AbstractTensor{T} where T <: Number
    data::DArray{T}
    ndim::Int8
    dims::Array{UInt32,1}
end
