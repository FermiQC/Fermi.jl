module ConfigurationInteraction

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
     ____                __  _                            _    _               
    / ___| ___   _ __   / _|(_)  __ _  _   _  _ __  __ _ | |_ (_)  ___   _ __  
   | |    / _ \ | '_ \ | |_ | | / _` || | | || '__|/ _` || __|| | / _ \ | '_ \ 
   | |___| (_) || | | ||  _|| || (_| || |_| || |  | (_| || |_ | || (_) || | | |
    \____|\___/ |_| |_||_|  |_| \__, | \__,_||_|   \__,_| \__||_| \___/ |_| |_|
    ___         _               |___/       _    _                             
   |_ _| _ __  | |_  ___  _ __  __ _   ___ | |_ (_)  ___   _ __                
    | | | '_ \ | __|/ _ \| '__|/ _` | / __|| __|| | / _ \ | '_ \
    | | | | | || |_|  __/| |  | (_| || (__ | |_ | || (_) || | | |
   |___||_| |_| \__|\___||_|   \__,_| \___| \__||_| \___/ |_| |_|
     
                                                                            
                          Module by G.J.R. Aroeira                           
================================================================================
"""
    @output "\n{}\n" banner
end

#defaults = Dict{Any,Any}(
#                :frozen => 0,
#                :active => -1,
#                :nroot => 1,
#                :min_matrix_elem => 1e-9
#)
#
#include("DetOperations.jl")
#include("MatrixElement.jl")
#include("FCI.jl")

end #module CI
