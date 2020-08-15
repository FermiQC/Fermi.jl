module Localize

using Lints
using Fermi.Orbitals: LocalOrbitals

function boys!(ints::IntegralHelper; thresh=1E-4) where LO 
    orbs = ints.orbs
    nao = length(orbs["*"][1].C)
end
end #module
