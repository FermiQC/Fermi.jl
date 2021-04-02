"""
    Fermi.HartreeFock

Module for running Hartree--Fock computations in Fermi.
"""
module HartreeFock
# Import Fermi basics
using Fermi
using Fermi.Options
using Fermi.Error
using Fermi: AbstractWavefunction
using Fermi.Geometry: Molecule
using Fermi.Integrals: IntegralHelper
using Fermi.Orbitals

function hf_header()
    output(repeat("=",80))
    output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:25}{:^28}{:25}|", "", "G.J.R Aroeira and M.M. Davis", "")
    output(repeat("=",80))
end

"""
    Fermi.HartreeFock.AbstractHFWavefunction

Abstract type common to all Hartree-Fock wave functions.

_struct tree:_

**AbstractHFWavefunction** <: AbstractWavefunction
"""
abstract type AbstractHFWavefunction <: AbstractWavefunction end

# Restricted Hartree--Fock
include("RHF/RHF.jl")

end #module

