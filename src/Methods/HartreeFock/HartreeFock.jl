"""
Module for running Hartree--Fock computations in Fermi.

"""
module HartreeFock

using Fermi
using Fermi.Integrals: ConventionalAOIntegrals
using Fermi.Geometry: Molecule
using Fermi.Output
using Lints
using LinearAlgebra

function print_header()
    @output repeat("=",80)*"\n"
    @output "|    {:<74}|\n" "Hartree Fock"
    @output "|        {:<70}|\n" "Module by M.M. Davis and G.J.R Aroeira"
    @output repeat("=",80)*"\n"
end

abstract type AbstractHFWavefunction <: Fermi.AbstractReferenceWavefunction end

include("RHF/RHF.jl")

end #module

