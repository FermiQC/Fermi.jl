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

abstract type AbstractOrbitals end

struct CanonicalOrbitals <: AbstractOrbitals
    orbs::Array{CanonicalOrbital}
end
struct NaturalOrbitals <: AbstractOrbitals
    orbs::Array{NaturalOrbitals}
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

end #module
