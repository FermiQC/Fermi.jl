"""
    Fermi.ConfigurationInteraction

Module for running ConfigurationInteraction computations in Fermi.
"""
module ConfigurationInteraction
# Import Fermi basics
using Fermi
using Fermi.Options
using Fermi.Integrals
using Fermi.Orbitals

function ci_header()
    banner = 
raw"""
================================================================================
//                       Configuration Interaction                            \\     
//                        Module by G.J.R. Aroeira                            \\       
================================================================================
"""
    output(banner)
end

"""
    Fermi.ConfigurationInteraction.AbstractCIWavefunction

Fermi abstract type common to all Configuration Interaction wavefunctions

_struct tree:_

**AbstractCIWavefunction** <: AbstractWavefunction
"""
abstract type AbstractCIWavefunction <: Fermi.AbstractWavefunction end

include("FCI/FCI.jl")

end #module CI
