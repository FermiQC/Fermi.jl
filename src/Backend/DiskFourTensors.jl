import Base.getindex
import Base.setindex!
import Base.+
import Base.eltype
struct DiskFourTensor{T}
    """
    Data structure for Rank 4 tensors
    """
    fname::String
    dname::String
    sz1::Int
    sz2::Int
    sz3::Int
    sz4::Int
    dtype::Type
end
function squeeze(A::AbstractArray)
    #singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
    return dropdims(A, dims = (findall(size(A) .== 1)...,))
end
function DiskFourTensor{T}(
    fname::String,
    sz1::Int,
    sz2::Int,
    sz3::Int,
    sz4::Int,
    mode::String = "r+",
) where T
    """
    Constructor for DiskFourTensor objects
    """
    file = h5open(fname, mode)
    if (mode == "w") | (mode == "w+")
        dataset = d_create(file, "data", datatype(dtype), dataspace(sz1))
    else
        dataset = file["data"]
    end
    close(file)
    DiskFourTensor(fname, "data", sz1, sz2, sz3, sz4, dtype)
end
# >>> overload getindex ( A[i,j,k,l] ) syntax

function getindex(
    dtens::DiskFourTensor,
    i1::Union{UnitRange{Int64},Int64,Colon},
    i2::Union{UnitRange{Int64},Int64,Colon},
    i3::Union{UnitRange{Int64},Int64,Colon},
    i4::Union{UnitRange{Int64},Int64,Colon},
)
    h5open(dtens.fname, "r") do fid
        squeeze(fid["$dtens.dname"][ranger(i1), ranger(i2), ranger(i3), ranger(i4)])
    end
end
# <<< 

# >>> overload setindex! ( A[i,j,k,l] = 2.0 )
function setindex!(
    dtens::DiskFourTensor,
    val,
    i1::UnitRange{Int64},
    i2::UnitRange{Int64},
    i3::UnitRange{Int64},
    i4::UnitRange{Int64},
)
    h5open(dtens.fname, "r+") do fid
        fid["$dtens.dname"][i1, i2, i3, i4] = val
    end
end
function setindex!(
    dtens::DiskFourTensor,
    val,
    i1::Union{UnitRange{Int64},Int64,Colon},
    i2::Union{UnitRange{Int64},Int64,Colon},
    i3::Union{UnitRange{Int64},Int64,Colon},
    i4::Union{UnitRange{Int64},Int64,Colon},
)
    h5open(dtens.fname, "r+") do fid
        fid["$dtens.dname"][ranger(i1), ranger(i2), ranger(i3), ranger(i4)] = val
    end
end
# <<< overload setindex!
# >>> overload eltype
function eltype(dtens::DiskFourTensor)
    return dtens.dtype
end
# <<< overload eltype

function blockfill!(dtens::DiskFourTensor, val)
    """
    Fill a DiskFourTensor with a single value.
    """
    A = zeros(Float64, dtens.sz1, dtens.sz2, dtens.sz3, dtens.sz4)
    A .= val
    h5write(dtens.fname, "$dtens.dname", A)
end
