"""
    DiskTensors

Disk based rank 1,2, and 4 tensors. Defines functions such as +, -, dot,... 
with buffered I/O. Uses HDF5 format for storing tensors. One tensor per file.

## structs
    DiskVector -> see ?Fermi.DiskTensors.DiskVector
    DiskMatrix -> see ?Fermi.DiskTensors.DiskMatrix
    DiskFourTensor -> see ?Fermi.DiskTensors.DiskFourTensor

## methods
    squeeze -> see ?Fermi.DiskTensors.squeeze
    eltype -> returns element type of a DiskFourTensor
    ranger -> converts input to a UnitRange, if possible
    blockfill! -> fills a DiskTensor with the value entered
"""
module DiskTensors

#include("DiskVectors.jl")
#include("DiskMatrices.jl")
include("DiskFourTensors.jl")

"""
	ranger
converts inp::Int64 -> UnitRange(inp:inp) for use in indexing DiskTensors
"""
function ranger(inp::Union{UnitRange{Int64},Int64,Colon})
    if typeof(inp) == Int64
        return UnitRange(inp:inp)
    else
        return inp
    end
end
export squeeze
export eltype
export ranger
#export DiskVector
export blockfill!
#export printdv
#export dvwrite!
#export dvread
#export dvdot
#export Base.getindex
#include("DiskVectorOps.jl")

#export DiskMatrix
#export dmread
#export dmwrite
#export dmdot

export DiskFourTensor
    #export dmelmult
end
