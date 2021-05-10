using LinearAlgebra
import Base: size, permutedims, getindex, setindex!, ndims, show, convert
import Base: iterate, length, similar, adjoint, eltype, +, -, *, /, ^, BroadcastStyle
import Strided: UnsafeStridedView

export FermiMDArray, FermiMDrand, FermiMDzeros, diagonalize

"""

    FermiMDArray{T,N}

Fermi array object held entirely in memory. Thin wrap around a standard Julia array, representing a dense array of type T and rank N.

_struct tree:_

**FermiMDArray** <: AbstractArray
"""
struct FermiMDArray{T,N} <: AbstractArray{T,N} 
    data::Array{T,N}
end

# Creates a FermiMDArray from a native Julia Array
function FermiMDArray(A::AbstractArray)
    return FermiMDArray(Array(A))
end

# If trying to create a FermiMDArray from another FermiMDArray, do nothing
function FermiMDArray(A::FermiMDArray{T,N}) where {T,N}
    return A
end

# If trying to create a FermiMDArray from a number, returns the number
function FermiMDArray(num::Number)
    return num
end

"""
    FermiMDzeros(x...)

Create a Julia Array as `zeros(x...)` wrapped in a FermiMDArray.
"""
function FermiMDzeros(x...)
    data = zeros(x...)
    return FermiMDArray(data)
end

"""
    FermiMDrand(x...)

Create a Julia Array as `rand(x...)` wrapped in a FermiMDArray.
"""
function FermiMDrand(x...)
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
    return size(A.data, i...)
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

"""
    Fermi.diagonalize(A::FermiMDArray; sortby=x->x)

Diagonalize a NxN FermiMDArray, returning eigenvalues and eigenvectors.

# Example
```
julia> A = FermiMDArray([3 6; 6 4])
julia> ϵ,ν = diagonalize(A)
julia> ϵ
Fermi Memory-held Dense Array - 2-element Array{Float64,1}:
 -2.5207972893961488
  9.520797289396146
```
"""
function diagonalize(A::FermiMDArray; sortby=x->x, hermitian=true)
    if hermitian
        vals, vecs = LinearAlgebra.eigen(Hermitian(A.data), sortby=sortby)
        return FermiMDArray(vals), FermiMDArray(vecs)
    else
        vals, vecs = LinearAlgebra.eigen(A.data, sortby=sortby)
        return FermiMDArray(vals), FermiMDArray(vecs)
    end
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

function UnsafeStridedView(A::FermiMDArray)
    UnsafeStridedView(A.data)
end

# Methods to allow Broadcasting
function Base.BroadcastStyle(::Type{<:FermiMDArray}) 
    Broadcast.ArrayStyle{FermiMDArray}()
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{FermiMDArray}}, ::Type{ElType}) where ElType
    FermiMDArray(similar(Array{ElType}, axes(bc)))
end

function Base.convert(::Type{T}, A::FermiMDArray) where T<:FermiMDArray{Float32}
    newdata = Base.convert(Array{Float32}, A.data)
    return FermiMDArray(newdata)
end

function Base.convert(::Type{T}, A::FermiMDArray) where T<:FermiMDArray{Float64}
    newdata = Base.convert(Array{Float64}, A.data)
    return FermiMDArray(newdata)
end

function show(io::IO, ::MIME"text/plain", A::FermiMDArray{T}) where T <: Number
    print("Fermi Memory-held Dense Array - ")
    display(A.data)
end