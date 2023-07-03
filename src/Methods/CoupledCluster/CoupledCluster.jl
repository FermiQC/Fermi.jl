"""
    Fermi.CoupledCluster

Module for running CoupledCluster computations in Fermi.
"""
module CoupledCluster
# Import Fermi basics
using Fermi
using Fermi.Options
using Fermi.Integrals
using Fermi.Orbitals

function cc_header()
    banner = 
raw"""
================================================================================
//    _____                   _          _   _____ _           _              \\ 
//   /  __ \                 | |        | | /  __ \ |         | |             \\   
//   | /  \/ ___  _   _ _ __ | | ___  __| | | /  \/ |_   _ ___| |_ ___ _ __   \\   
//   | |    / _ \| | | | '_ \| |/ _ \/ _` | | |   | | | | / __| __/ _ \ '__|  \\   
//   | \__/\ (_) | |_| | |_) | |  __/ (_| | | \__/\ | |_| \__ \ ||  __/ |     \\  
//    \____/\___/ \__,_| .__/|_|\___|\__,_|  \____/_|\__,_|___/\__\___|_|     \\  
//                    | |                                                     \\   
//                    |_|                                                     \\   
//                                                                            \\     
//                 Module by G.J.R. Aroeira and M. M. Davis                   \\       
================================================================================
"""
    output(banner)
end

"""
    Fermi.CoupledCluster.AbstractCCWavefunction

Fermi abstract type common to all Coupled Cluster wavefunctions

_struct tree:_

**AbstractCCWavefunction** <: AbstractWavefunction
"""
abstract type AbstractCCWavefunction <: Fermi.AbstractWavefunction end

#include("RCCD/RCCD.jl")
include("RCCSD/RCCSD.jl")
include("PerturbativeTriples/PerturbativeTriples.jl")
include("ExternallyCorrected/ecCCSD/ecRCCSD.jl")

end #module CC
