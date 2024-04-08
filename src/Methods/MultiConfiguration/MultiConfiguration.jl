"""
    Fermi.MultiConfiguration

Module for running Multiconfigurational computations in Fermi.

# Methods

    > Fermi.MultiConfiguration.MCSCF
    > Fermi.MultiConfiguration.MRPT2
"""
module MultiConfiguration
# Import Fermi basics
using LinearAlgebra: hermitian, vcat
using Fermi
using Fermi.Options
using Fermi.Integrals
using Fermi.Orbitals

function mc_header()
    output(repeat("=",80))
    output("|{:22}{:^32}{:22}|", "", "Multiconfigurational Calculation", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:36}{:^8}{:36}|", "", "S. Jeong", "")
    output(repeat("=",80))
end

# function umc_header()
#     output(repeat("=",80))
#     output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
#     output("|{:34}{:^9}{:34}|", "", "Module  by","")
#     output("|{:17}{:^44}{:17}|", "", "G.J.R Aroeira, M.M. Davis, and S.M. Goodlett", "")
#     output(repeat("=",80))
# end

"""
    Fermi.HartreeFock.AbstractHFWavefunction

Abstract type common to all Hartree-Fock wave functions.

Struct tree

**AbstractHFWavefunction** <: AbstractWavefunction
"""
abstract type AbstractMCWavefunction <: Fermi.AbstractWavefunction end

# Different MCSCF methods are included here:
# Multiconfigurational SCF calculations
include("MCSCF/MCSCF.jl")

# Multireference Perturbation Theory calculations
# include("MRPT/MRPT.jl")

# Multireference Perturbation Theory calculations
# include("MRCI/MRCI.jl")

end #module
