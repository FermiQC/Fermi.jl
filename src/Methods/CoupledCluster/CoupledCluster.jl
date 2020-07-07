"""
Module for running CC computations in Julia.

Implemented --> RCCD, RCCSD, DF-RCCD
"""
module CoupledCluster

using Fermi.Wavefunction
using Fermi.Transformation
using Fermi.Output
using Printf
using Base.Threads
using TensorOperations
using LinearAlgebra
using Dates

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

#defaults = Dict(
#                :cc_max_iter => 50,
#                :cc_max_rms => 10^-10,
#                :cc_e_conv => 10^-10,
#                :diis => false,
#                :do_pT => false,
#                :fcn => 0
#               )
#
#include("PerturbativeTriples.jl")
#include("RCCD.jl")
#include("DF-RCCD.jl")
#include("ROCCD.jl")
#include("RCCSD.jl")
#include("UCCSD.jl")
#include("mRCCD.jl")
#include("AutoRCCSD.jl")
#include("ecRCCSD.jl")

end #module CC
