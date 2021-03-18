"""
    Fermi.MøllerPlesset

Module for running Møller--Plesset perturbation theory computations.
"""

module MollerPlesset
# Import Fermi basics
using Fermi
using Fermi.Options
using Fermi.Error
using Fermi: AbstractWavefunction
using Fermi.Geometry: Molecule
using Fermi.Integrals

function mp_header()
    output("\n================================================================================")
    output("|   Møller-Plesset Perturbation Theory                                         |")
    output("|       Module written by G.J.R. Aroeira and M.M. Davis                        |")
    output("================================================================================") 
end

abstract type AbstractMPWavefunction <: Fermi.AbstractWavefunction end

include("RMP2/RMP2.jl")
end #module
