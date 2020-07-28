module ConfigurationInteraction

using Fermi
using Fermi.Geometry: Molecule
using Fermi.Integrals: ConventionalAOIntegrals
using Fermi.Output
using TensorOperations
using LinearAlgebra

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
     
                                                                            
             Module by G.J.R. Aroeira and M.M. Davis                        
================================================================================
"""
    @output "\n{}\n" banner
end

abstract type AbstractCIWavefunction <: Fermi.AbstractCorrelatedWavefunction end

# Struct symbolizing the type of implementation for different CI methods
abstract type CIAlgorithm end
struct SparseHamiltonian <: CIAlgorithm end
struct ACI <: CIAlgorithm end


include("CASCI/CASCI.jl")

end #module CI
