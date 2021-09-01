"""
    Fermi.HartreeFock

Module for running Hartree--Fock computations in Fermi.

# Methods

    > Fermi.HartreeFock.RHF
    > Fermi.HartreeFock.UHF
"""
module HartreeFock
# Import Fermi basics
using LinearAlgebra: hermitian, vcat
using Fermi
using Fermi.Options
using Fermi.Integrals
using Fermi.Orbitals

function hf_header()
    output(repeat("=",80))
    output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:25}{:^28}{:25}|", "", "G.J.R Aroeira and M.M. Davis", "")
    output(repeat("=",80))
end

function uhf_header()
    output(repeat("=",80))
    output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:17}{:^44}{:17}|", "", "G.J.R Aroeira, M.M. Davis, and S.M. Goodlett", "")
    output(repeat("=",80))
end

"""
    Fermi.HartreeFock.AbstractHFWavefunction

Abstract type common to all Hartree-Fock wave functions.

Struct tree

**AbstractHFWavefunction** <: AbstractWavefunction
"""
abstract type AbstractHFWavefunction <: Fermi.AbstractWavefunction end

"""
    Fermi.HartreeFock.get_scf_alg()

Returns a singleton type corresponding to a RHF implementation based on the options.
"""

function get_scf_alg()
    ref = Options.get("reference")
    if ref == "rhf"
        implemented = [RHFa()]
    elseif ref == "uhf"
        implemented = [UHFa()]
    else
        throw(FermiException("$ref is not a supported reference energy"))
    end

    N = Options.get("scf_alg")
    try
        return implemented[N]
    catch BoundsError
        throw(FermiException("implementation number $N not available for $ref."))
    end
end

# Different Hartree-Fock methods are included here:

# Restricted Hartree--Fock
include("RHF/RHF.jl")

# UHF
include("UHF/UHF.jl")

end #module
