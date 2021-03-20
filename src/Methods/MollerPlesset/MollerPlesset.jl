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
using Fermi.Integrals: IntegralHelper

function mp_header()
    output(repeat("=",80))
    output("|{:22}{:^12}{:22}|", "Møller-Plesset Perturbation Theory", "", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:25}{:^28}{:25}|", "", "G.J.R Aroeira and M.M. Davis", "")
    output(repeat("=",80))
end

"""
    Fermi.MollerPlesset.AbstractMPWavefunction

Abstract type common to all Møller-Plesset wave functions.

_struct tree:_

**AbstractMPWavefunction** <: AbstractWavefunction
"""
abstract type AbstractMPWavefunction <: Fermi.AbstractWavefunction end

# Restricted MP2
include("RMP2/RMP2.jl")

end #module
