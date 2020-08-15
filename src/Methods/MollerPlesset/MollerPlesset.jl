"""
    Fermi.MollerPlesset

module for running MP2 energies on restricted and unrestricted HF references.
## methods
    do_rmp2 -> see docstring ?Fermi.MollerPlesset.do_rmp2
    do_ump2 -> see docstring ?Fermi.MollerPlesset.do_ump2
"""
module MollerPlesset

using Fermi
#using Fermi.Direct
using Fermi.Output
using TensorOperations

export RMP2
export RMP3

function print_header()
    @output "================================================================================\n" 
    @output "|   Moller-Plesset Perturbation Theory                                         |\n"
    @output "|       module written by M.M. Davis                                           |\n"
    @output "================================================================================\n" 
end

abstract type MP2Algorithm end

struct Conventional <: MP2Algorithm end
struct DF           <: MP2Algorithm end
struct Direct       <: MP2Algorithm end

abstract type AbstractMPWavefunction <: Fermi.AbstractCorrelatedWavefunction end

include("RMP2.jl")
include("RMP3.jl")

end #module
