module DIIS
using LinearAlgebra
import Base.push!
import Base.length

struct DIISManager{T1<:AbstractFloat,
                  T2 <: AbstractFloat }# where T <: AbstractFloat
    vecs::Array{Array{T1},1}
    errs::Array{Array{T2},1}
    max_vec::Int64
end

struct CROPManager{T<:AbstractFloat}
    Wopt::Array{Array{T}}
    Topt::Array{Array{T}}
    Waux::Array{Array{T}}
    Taux::Array{Array{T}}
end



function DIISManager{T1,T2}(;size=6) where { T1 <: AbstractFloat,
                                             T2 <: AbstractFloat }
    vecs = Array{Array{T1}}(undef,0)
    errs = Array{Array{T2}}(undef,0)
    DIISManager{T1,T2}(vecs,errs,size)
end

function CROPManager{T}(;size=6) where T
    Wopt = Array{Array{T}}(undef,0)
    Topt = Array{Array{T}}(undef,0)
    Waux = Array{Array{T}}(undef,0)
    Taux = Array{Array{T}}(undef,0)
    CROPManager{T}(Wopt,Topt,Waux,Taux)
end

function extrap(M::CROPManager{T}) where T
end
    
function length(M::DIISManager{T1,T2}) where { T1 <: AbstractFloat,
                                               T2 <: AbstractFloat }
    length(M.vecs)
end

function push!(M::DIISManager{T1,T2}, V::Array, E::Array) where { T1 <: AbstractFloat,
                                                                  T2 <: AbstractFloat }
    if length(M)+1 > M.max_vec
        norms = norm.(M.errs)
        idx = findmax(norms)[2]
        deleteat!(M.vecs,idx)
        deleteat!(M.errs,idx)
    end
    push!(M.vecs,convert(Array{T1},deepcopy(V)))
    push!(M.errs,convert(Array{T2},deepcopy(E)))
end

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
    B[1:E-1,1:E-1] ./= maximum(abs.(B[1:E-1,1:E-1]))
    resid = zeros(T1,diis_size+1)
    resid[end] = 1
    LAPACK.gesv!(B,resid)
    ci = resid
    out = zeros(T1,size(M.vecs[1]))
    outw = zeros(T1,size(M.vecs[1]))
    for num in 1:diis_size
        out += ci[num]*M.vecs[num]
        outw += ci[num]*M.errs[num]
    end
    add_res ? out += outw : nothing
    out
end
end #module
