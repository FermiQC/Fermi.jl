module DIIS
using LinearAlgebra
import Base.push!
import Base.length

"""
    DIISManager{T1<:AbstractFloat,
                T2<:AbstractFloat}

# Fields

    vecs::Array{Array{T1},1} vectors to extrapolate from.
    errs::Array{Array{T2},1} error vectors to build B matrix from.
    max_vec::Int64           max number of vectors to hold
"""
struct DIISManager{T1<:AbstractFloat,T2 <: AbstractFloat}
    vecs::Array{AbstractArray{T1},1}
    errs::Array{AbstractArray{T2},1}
    max_vec::Int64
end

struct CROPManager{T<:AbstractFloat}
    Wopt::Array{AbstractArray{T}}
    Topt::Array{AbstractArray{T}}
    Waux::Array{AbstractArray{T}}
    Taux::Array{AbstractArray{T}}
end

function DIISManager{T1,T2}(;size=6) where { T1 <: AbstractFloat,
                                             T2 <: AbstractFloat }
    vecs = Array{AbstractArray{T1}}(undef,0)
    errs = Array{AbstractArray{T2}}(undef,0)
    DIISManager{T1,T2}(vecs,errs,size)
end

function CROPManager{T}(;size=6) where T
    Wopt = Array{Array{T}}(undef,0)
    Topt = Array{Array{T}}(undef,0)
    Waux = Array{Array{T}}(undef,0)
    Taux = Array{Array{T}}(undef,0)
    CROPManager{T}(Wopt,Topt,Waux,Taux)
end

function length(M::DIISManager) 
    length(M.vecs)
end

function push!(M::DIISManager{T1,T2}, V::AbstractArray, E::AbstractArray) where { T1 <: AbstractFloat,
                                                                  T2 <: AbstractFloat }
    if length(M)+1 > M.max_vec
        norms = norm.(M.errs)
        idx = findmax(norms)[2]
        deleteat!(M.vecs,idx)
        deleteat!(M.errs,idx)
    end
    push!(M.vecs,deepcopy(V))
    push!(M.errs,deepcopy(E))
end

"""
    extrapolate(M::DIISManager{T1,T2}; add_res=false) where { T1 <: AbstractFloat,
                                                              T2 <: AbstractFloat }

Takes current state of M and produces an optimal (within subspace) trial vector.

# kwargs
    add_res=false      Do add estimated residual to trial vector? 
"""
function extrapolate(M::DIISManager{T1,T2};add_res=false) where { T1 <: AbstractFloat,
                                                    T2 <: AbstractFloat }
    diis_size = length(M)
    B = ones(T1,diis_size+1,diis_size+1)*1
    B[end,end] = 0
    for (n1, e1) in enumerate(M.errs[1:end])
        for (n2, e2) in enumerate(M.errs[1:end])
            B[n1,n2] = sum(e1 .* e2) #(e1,e2)
        end
    end 
    E = size(B,1)
    resid = zeros(T1,diis_size+1)
    resid[end] = 1
    ci = svd(B)\resid
    out = zero(M.vecs[1])
    add_res ? outw = zeros(T1,size(M.vecs[1])) : nothing
    for num in 1:diis_size
        out += ci[num]*M.vecs[num]
        add_res ? outw += ci[num]*M.errs[num] : nothing
    end
    add_res ? out += outw : nothing
    out
end
end #module
