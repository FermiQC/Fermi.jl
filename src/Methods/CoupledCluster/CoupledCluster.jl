"""
Module for running CoupledCluster computations in Fermi.

"""
module CoupledCluster

using Fermi
using Fermi.Integrals: IntegralHelper
using Fermi.HartreeFock: RHF
using Fermi.Output
using TensorOperations
using LinearAlgebra

function print_header()
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
//                 Module by M.M. Davis and G.J.R. Aroeira                    \\       
================================================================================
"""
    @output "\n{}\n" banner
end

"""
    Fermi.CoupledCluster.AbstractCCWavefunction

Fermi abstract type common to all Coupled Cluster wavefunctions

_struct tree:_

**AbstractCCWavefunction** <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
abstract type AbstractCCWavefunction <: Fermi.AbstractCorrelatedWavefunction end

# Structures symbolizing the type of implementation for each CC method
abstract type CCAlgorithm end
struct DPD   <: CCAlgorithm end
struct CTF   <: CCAlgorithm end
struct DFCTF <: CCAlgorithm end

function select_precision(A::String)
    implemented = Dict{Any,Any}("single" => Float32,
                                "double" => Float64)
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid precision: $(A)"))
    end
end

function select_algorithm(A::String)
    implemented = Dict{Any,Any}("DPD" => DPD(),
                                "CTF" => CTF(),
                                "DFCTF" => DFCTF())
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid CC algorithm: $(A)"))
    end
end

#include("RCCD/RCCD.jl")
include("RCCSD/RCCSD.jl")
include("pT/PerturbativeTriples.jl")
include("ecRCCSD/ecRCCSD.jl")
include("ecpT/ecPerturbativeTriples.jl")
include("BCC/BD.jl")

end #module CC
