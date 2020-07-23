module DIIS
using LinearAlgebra
import Base.push!
import Base.length

struct DIISManager{T<:AbstractFloat}# where T <: AbstractFloat
    vecs::Array{Array{T},1}
    errs::Array{Array{T},1}
    max_vec::Int64
end

function DIISManager{T}(;size=6) where T <: AbstractFloat
    vecs = Array{Array{T}}(undef,0)
    errs = Array{Array{T}}(undef,0)
    DIISManager{T}(vecs,errs,size)
end
    
function length(M::DIISManager{T}) where T <: AbstractFloat
    length(M.vecs)
end

function push!(M::DIISManager{T}, V::Array, E::Array) where T <: AbstractFloat
    push!(M.vecs,deepcopy(V))
    push!(M.errs,deepcopy(E))
    #push!(M.errs,reshape(deepcopy(E),prod(size(V))))
    if length(M) > M.max_vec
        norms = norm.(M.errs)
        idx = findmax(norms)[2]
        deleteat!(M.vecs,idx)
        deleteat!(M.errs,idx)
    end
end

function extrapolate(M::DIISManager{T}) where T <: AbstractFloat
    diis_size = length(M)
    B = ones(T,diis_size+1,diis_size+1)*1
    B[end,end] = 0
    for (n1, e1) in enumerate(M.errs[1:end])
        for (n2, e2) in enumerate(M.errs[1:end])
            B[n1,n2] = sum(e1 .* e2) #(e1,e2)
        end
    end 
    E = size(B,1)
    B[1:E-1,1:E-1] ./= maximum(abs.(B[1:E-1,1:E-1]))
    resid = zeros(T,diis_size+1)
    resid[end] = 1
    LAPACK.gesv!(B,resid)
    ci = resid
    out = zeros(T,size(M.vecs[1]))
    for num in 1:diis_size
        out .+= ci[num]*M.vecs[num]
    end
    out
end
end #module
