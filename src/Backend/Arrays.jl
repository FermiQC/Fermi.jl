using LinearAlgebra
import Base: size, permutedims, getindex, setindex!, ndims, show, iterate, length, similar, adjoint, eltype, +, -, *, /, ^, BroadcastStyle

export FermiMDArray, FermiMDRand, FermiMDZeros, diagonalize

"""

    FermiMDArray{T,N}

Fermi array concrete object held entirely in memory. Thin wrapper around a standard Julia array, representing a dense array of type T and rank N.

_struct tree:_

**FermiMDArray** <: AbstractArray
"""
struct FermiMDArray{T,N} <: AbstractArray{T,N}
    data::AbstractArray{T,N}
end

function FermiMDArray(data::AbstractArray{T,N}) where {T,N}
    eltype(data) <: Number || throw(TypeError(:FermiMDArray, "data type: must be numerical", Number, eltype(data)))
    FermiMDArray{T,N}(data)
end

function FermiMDArray(num::Number)
    return num
end

# Quick constructors for Fermi Arrays
function FermiMDZeros(x...)
    data = zeros(x...)
    return FermiMDArray(data)
end

function FermiMDRand(x...)
    data = rand(x...)
    return FermiMDArray(data)
end

function similar(A::FermiMDArray)
    data = similar(A.data)
    return FermiMDArray(data)
end

function similar(A::FermiMDArray, dims::Dims)
    data = similar(A.data, dims)
    return FermiMDArray(data)
end

# Basic methods for AbstractArrays in Julia
function size(A::FermiMDArray, i...)
    size(A.data, i...)
end

function getindex(A::FermiMDArray, I...)
    return FermiMDArray(A.data[I...])
end

function setindex!(A::FermiMDArray, val, I...)
    A.data[I...] = val
end

function ndims(A::FermiMDArray)
    return ndims(A.data)
end

function length(A::FermiMDArray)
    return length(A.data)
end

function eltype(A::FermiMDArray)
    return eltype(A.data)
end

# Special methods for LinearAlgebra
function permutedims(A::FermiMDArray,tup)
    FermiMDArray(permutedims(A.data,tup))
end

function permutedims!(A::FermiMDArray,tup)
    A.data .= permutedims(A.data,tup)
end

function diagonalize(A::FermiMDArray; sortby=x->x)
    vals, vecs = LinearAlgebra.eigen(A.data, sortby=sortby)
    return FermiMDArray(vals), FermiMDArray(vecs)
end

function LinearAlgebra.Hermitian(A::FermiMDArray)
    return FermiMDArray(LinearAlgebra.Hermitian(A.data))
end

function LinearAlgebra.diag(A::FermiMDArray)
    return FermiMDArray(LinearAlgebra.diag(A.data))
end

function adjoint(A::FermiMDArray)
    return FermiMDArray(adjoint(A.data))
end

# Basic mathematical methods

function Base.:+(A::FermiMDArray, B::FermiMDArray)
    return FermiMDArray(A.data + B.data)
end

function Base.:+(A::FermiMDArray, B::AbstractArray)
    return FermiMDArray(A.data + B)
end

function Base.:+(A::AbstractArray, B::FermiMDArray)
    return FermiMDArray(A + B.data)
end

function Base.:-(A::FermiMDArray, B::FermiMDArray)
    return FermiMDArray(A.data - B.data)
end

function Base.:-(A::FermiMDArray, B::AbstractArray)
    return FermiMDArray(A.data - B)
end

function Base.:-(A::AbstractArray, B::FermiMDArray)
    return FermiMDArray(A - B.data)
end

function Base.:*(A::FermiMDArray, B::Number)
    return FermiMDArray(B*A.data)
end

function Base.:*(A::Number, B::FermiMDArray)
    return FermiMDArray(A*B.data)
end

function Base.:*(A::FermiMDArray, B::FermiMDArray)
    return FermiMDArray(A.data*B.data)
end

function Base.:^(A::FermiMDArray, B::Integer)
    return FermiMDArray(A.data^B)
end

function Base.:^(A::FermiMDArray, B::AbstractFloat)
    return FermiMDArray(A.data^B)
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", A::FermiMDArray{T}) where T <: Number
    print("Fermi Memory-held Dense Array - ")
    display(A.data)
end

# Methods to allow Broadcasting
function Base.BroadcastStyle(::Type{<:FermiMDArray}) 
    Broadcast.ArrayStyle{FermiMDArray}()
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{FermiMDArray}}, ::Type{ElType}) where ElType
    FermiMDArray(similar(Array{ElType}, axes(bc)))
end
