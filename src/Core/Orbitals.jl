module Orbitals

import Base.getindex
import Base.lastindex

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

function getindex(orb::O,i) where O <: AbstractOrbitals
    return orb.orbs[i]
end
function getindex(orb::O,R::UnitRange) where O <: AbstractOrbitals
    return O([orb.orbs[i] for i in R])
end

function lastindex(orb::O) where O <: AbstractOrbitals
    return length(orb.orbs)
end

function rotate!(orbs::O,U::Array{Float64,2}) where O <: AbstractOrbitals
    C = hcat([orb.C for orb in orbs.orbs]...)
    C = C*U
    for i in eachindex(orbs.orbs)
        orbs[i].C[:] .= C[:,i]
    end
end

#mutable struct OrbDict
#    all::O where O <: AbstractOrbitals
#    ndocc::Int
#    nvir::Int
#end
#
#function getindex(OD::OrbDict,entry::String)
#    if entry == "*"
#        return OD.all
#    elseif entry = "O" || entry = "o"
#        return CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in 1:ndocc])
#    elseif entry = "V" || entry = "v"
#        return CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in ndocc+1:ndocc+nvir])
#    end
#end


end #module
