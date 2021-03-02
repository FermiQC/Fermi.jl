"""
    FermiArray

Abstract type for arrays used within Fermi

_struct tree:_

**FermiArray**
"""
abstract type FermiArray end

import Base: size, permutedims, getindex, ndims, show, iterate, length, similar, eltype

export FermiMDArray, FermiMDRand, FermiMDZeros

"""

    FermiMDArray{T}

Fermi array concrete object held entirely in memory. Thin wrapper around a standard Julia array, representing a dense array of type T, where T <: Number.

_struct tree:_

**FermiMDArray** <: FermiArray
"""
struct FermiMDArray{T} <: FermiArray where T <: Number
    data::Array{T}
end

# Type check
function FermiMDArray(data::Array)
    eltype(data) <: Number || throw(TypeError(:FermiMDArray, "data type: must be numerical", Number, eltype(data)))
    FermiMDArray{eltype(data)}(data)
end

function FermiMDZeros(x...)
    data = zeros(x...)
    return FermiMDArray(data)
end

function FermiMDRand(x...)
    data = rand(x...)
    return FermiMDArray(data)
end

# Methods for typical operations with arrays in Julia
function size(A::FermiMDArray, i...)
    size(A.data, i...)
end

function permutedims(A::FermiMDArray,tup)
    FermiMDArray(permutedims(A.data,tup))
end

function permutedims!(A::FermiMDArray,tup)
    A.data .= permutedims(A.data,tup)
end

function getindex(A::FermiMDArray,I...)
    return A.data[I...]
end

function ndims(A::FermiMDArray)
    return ndims(A.data)
end

function length(A::FermiMDArray)
    return length(A.data)
end

function iterate(A::FermiMDArray, i...)
    return iterate(A.data, i...)
end

function similar(A::FermiMDArray, i...)
    data = similar(A.data, i...)
    return FermiMDArray(data)
end

function eltype(A::FermiMDArray)
    return eltype(A.data)
end

function show(io::IO, ::MIME"text/plain", A::FermiMDArray{T}) where T <: Number
    print("Fermi Memory-held Dense Array - ")
    display(A.data)
end
