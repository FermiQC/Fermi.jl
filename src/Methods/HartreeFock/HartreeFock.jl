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

export RHF

function hf_header()
    output(repeat("=",80))
    output("|    {:<74}|\n", "Hartree Fock", ending="")
    output("|        {:<70}|\n", "Module by M.M. Davis and G.J.R Aroeira", ending="")
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

