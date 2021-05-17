module DIIS
using LinearAlgebra
import Base.push!
import Base.length

export DIISManager

"""
    DIISManager{T1<:AbstractFloat,
                T2<:AbstractFloat}

Structure holding arrays necessary for extrapolation with the `DIIS` method.

# Fields

    vecs                    Vectors to extrapolate from
    errs                    Error vectors to build B matrix from
    max_vec::Int64          Max number of vectors to hold
"""
struct DIISManager{T1<:AbstractFloat,T2 <: AbstractFloat}
    vecs::Array{AbstractArray{T1},1}
    errs::Array{AbstractArray{T2},1}
    max_vec::Int64
end

function DIISManager{T1,T2}(;size=6) where { T1 <: AbstractFloat,
                                             T2 <: AbstractFloat }
    vecs = Array{AbstractArray{T1}}(undef,0)
    errs = Array{AbstractArray{T2}}(undef,0)
    DIISManager{T1,T2}(vecs,errs,size)
end

function length(M::DIISManager) 
    length(M.vecs)
end

function push!(M::DIISManager{T1,T2}, V::AbstractArray, E::AbstractArray) where { T1 <: AbstractFloat,
                                                                  T2 <: AbstractFloat }
    # If the number of vectors stored has reach its maximum
    # the vector associated with the biggest error is replaced with the new one
    if length(M)+1 > M.max_vec
        norms = norm.(M.errs)
        _, idx = findmax(norms)
        deleteat!(M.vecs,idx)
        deleteat!(M.errs,idx)
    end
    push!(M.vecs,deepcopy(V))
    push!(M.errs,deepcopy(E))
end

"""
    Fermi.DIIS.extrapolate(M::DIISManager{T1,T2}; add_res=false) where { T1 <: AbstractFloat,
                                                              T2 <: AbstractFloat }

Produces a new guess vector using *direct inversion in the iterative subspace* given 
the information in `M` which is a `DIISManager`.
"""
function extrapolate(M::DIISManager{T1,T2}; add_res=false) where {T1 <: AbstractFloat, T2 <: AbstractFloat}

    if length(M) == 1
        throw(DIISError(" cannot extrapolate from one vector"))
    end

    # Solves the equation for the new vector
    N = length(M) + 1

    # Create the B matrix
    B = ones(T1, N, N)
    B[end,end] = zero(T1)
    for (n1, e1) in enumerate(M.errs[1:end])
        for (n2, e2) in enumerate(M.errs[1:end])
            # Bij = (ei⋅ej)
            B[n1,n2] = sum(e1 .* e2) 
        end
    end 

    # Create the residual vector
    resid = zeros(T1, N)
    resid[end] = one(T1) 

    # Get coefficients
    ci = B\resid
    out = zero(M.vecs[1])

    # Compute the new vector as Pnew = ∑ ci * Pi
    add_res ? outw = zeros(T1,size(M.vecs[1])) : nothing
    out = sum(ci[1:end-1] .* M.vecs)    # Need to pop the last element as it is λ not ci

    # Add this extra stuff if requested
    if add_res
        out += sum(ci[1:end-1] .* M.errs)
    end

    return out
end
end #module
