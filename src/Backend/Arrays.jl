"""
    FermiArray

Abstract type for arrays used within Fermi

_struct tree:_

**FermiArray**
"""
abstract type FermiArray end

import Base.size
import Base.permutedims
import Base.getindex
import Base: ndims

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

function size(A::FermiMDArray)
    size(A.data)
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
