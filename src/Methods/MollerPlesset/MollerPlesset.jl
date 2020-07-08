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
using Fermi.IntegralTransformation
using Fermi.Output
using TensorOperations

export do_rmp2
export do_ump2
export do_direct_rmp2
export do_df_rmp2
export RMP2

abstract type AbstractMPWavefunction <: Fermi.AbstractCorrelatedWavefunction end

function print_header()
    @output "================================================================================\n" 
    @output "|   Moller-Plesset Perturbation Theory                                         |\n"
    @output "|       module written by M.M. Davis                                           |\n"
    @output "================================================================================\n" 
end

include("RMP2.jl")
#include("UMP2.jl")
#include("DirectRMP2.jl")
#include("DF-RMP2.jl")

end #module
