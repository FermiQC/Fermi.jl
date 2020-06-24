module SpecialArrays

using Base: getindex

struct T3αβα
    o::Int64
    v::Int64
    A::Array{Array{Array{Float64,2},1},1}
end


