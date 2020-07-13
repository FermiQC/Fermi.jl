"""
Module for running Hartree--Fock computations in Fermi.

"""
module HartreeFock

using Fermi
using Fermi.Output
using Lints

function print_header()
    @output repeat("=",80)*"\n"
    @output "|    {:<74}|\n" "Hartree Fock"
    @output "|        {:<70}|\n" "Module by M.M. Davis and G.J.R Aroeira"
    @output repeat("=",80)*"\n"
end

abstract type AbstractHFWavefunction <: Fermi.AbstractReferenceWavefunction end

include("RHF/RHF.jl")

end #module

