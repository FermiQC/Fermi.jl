"""
    Fermi.HartreeFock

Module for running Hartree--Fock computations in Fermi.
"""
module HartreeFock

using Fermi.Output

function print_header()
    @output repeat("=",80)*"\n"
    @output "|    {:<74}|\n" "Hartree Fock"
    @output "|        {:<70}|\n" "Module by M.M. Davis and G.J.R Aroeira"
    @output repeat("=",80)*"\n"
end

"""
    Fermi.HartreeFock.AbstractHFWavefunction

Abstract type common to all Hartree-Fock wave functions.

_struct tree:_

**AbstractHFWavefunction** <: AbstractWavefunction
"""
abstract type AbstractHFWavefunction <: Fermi.AbstractWavefunction end

# Restricted Hartree--Fock
include("RHF/RHF.jl")

end #module

