module Orbitals

import Base.getindex
import Base.lastindex
import Base.setindex!
import Base.length

abstract type AbstractOrbital end

struct CanonicalOrbital <: AbstractOrbital
    C::Array{Float64,1}
end # module
struct NaturalOrbital <: AbstractOrbital
    C::Array{Float64,1}
    ω::Float64 # occupation number
end
struct LocalOrbital <: AbstractOrbital
    C::Array{Float64,1}
end
struct BruecknerOrbital <: AbstractOrbital
    C::Array{Float64,1}
end

abstract type AbstractOrbitals end

struct CanonicalOrbitals <: AbstractOrbitals
    orbs::Array{CanonicalOrbital}
end
struct NaturalOrbitals <: AbstractOrbitals
    orbs::Array{NaturalOrbitals}
end
struct LocalOrbitals <: AbstractOrbitals
    orbs::Array{LocalOrbitals}
end
struct BruecknerOrbitals <: AbstractOrbitals
    orbs::Array{BruecknerOrbital}
end

function CanonicalOrbitals(C::Array{Float64,2})
    return CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in 1:size(C,2)])
end

function toarray(orbs::O) where O <: AbstractOrbitals
    C = hcat([orb.C for orb in orbs.orbs]...)
end

function Base.getindex(orb::O,i) where O <: AbstractOrbitals
    return orb.orbs[i]
end
function Base.getindex(orb::O,R::UnitRange) where O <: AbstractOrbitals
    return O([orb.orbs[i] for i in R])
end

function Base.lastindex(orb::O) where O <: AbstractOrbitals
    return length(orb.orbs)
end

function Base.length(orb::O) where O <: AbstractOrbitals
    return length(orb.orbs)
end


mutable struct OrbDict
    current::String
    ndocc::Int
    nvir::Int
    frozencore::Int
    frozenvir::Int
    cache::Dict{String,O} where O <: AbstractOrbitals
end

function matchphase(orbs,oldC)
    for oi in eachindex(orbs.cache[orbs.current].orbs)
        o = orbs.cache[orbs.current].orbs[oi]
        osorted = reverse(sort(o.C,by=abs))
        csorted = reverse(sort(oldC[:,oi],by=abs))
        if sign(osorted[1]) != sign(csorted[1])
            o.C .*= -1
        end
    end
end
function topositive!(orbs::OrbDict)
    for o in orbs.cache[orbs.current].orbs
        sorted = reverse(sort(o.C,by=abs))
        if sorted[1] >= 0
            o.C .*= 1
        else
            o.C .*= -1
        end
    end
end
function rotate!(orbs::OrbDict,U::Array{Float64,2};fc=0,fv=0) where O <: AbstractOrbitals
    C = hcat([orb.C for orb in orbs["OV"].orbs]...)
    C = C*U
    for i=1:length(orbs["OV"])
        _i = i+fc
        orbs.cache[orbs.current].orbs[_i].C[:] .= C[:,i]
    end
end

function OrbDict()
    OrbDict("canonical",0,0,0,0,Dict{String,AbstractOrbitals}())
end

function activate!(OD::OrbDict,entry::String)
    if !(entry in keys(OD.cache))
        error("OrbDict tried to activate a non-existent basis!")
    else
        OD.current = entry
    end
end

function Base.getindex(OD::OrbDict,entry::String)
    gather = occursin("[",entry)
    s = strip(entry,['[',']'])
    all = OD.cache[OD.current]
    type = typeof(all)
    fc = OD.frozencore
    fv = OD.frozenvir
    lookup = Dict{Char,type}( 
                             'F' => all[1:OD.ndocc],
                             'O' => all[fc+1:OD.ndocc],
                             'U' => all[OD.ndocc+1:OD.ndocc+OD.nvir],
                             'V' => all[OD.ndocc+1:OD.ndocc+OD.nvir-fv]
                            )
    out = []
    for i in eachindex(s)
        _s = lookup[s[i]]
        if length(_s) != 0
            push!(out,lookup[s[i]])
        end
    end
    out = type(hcat([toarray(o) for o in out]...))
    if gather
        return toarray(out)
    else
        return out
    end
end

function Base.setindex!(OD::OrbDict,val,entry::String)
    OD.cache[entry] = val
end

end #module
