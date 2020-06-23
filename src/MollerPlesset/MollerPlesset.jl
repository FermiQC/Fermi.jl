"""
    Fermi.MollerPlesset

module for running MP2 energies on restricted and unrestricted HF references.
## methods
    do_rmp2 -> see docstring ?Fermi.MollerPlesset.do_rmp2
    do_ump2 -> see docstring ?Fermi.MollerPlesset.do_ump2
"""
module MollerPlesset

using Fermi.Wavefunction
using Fermi.DiskTensors
#using Fermi.Direct
using Fermi.IntegralTransformation
using Fermi.DF
using Fermi.Output
using Fermi
using TensorOperations

export do_rmp2
export do_ump2
export do_direct_rmp2
export do_df_rmp2

function print_header()
    @output "================================================================================\n" 
    @output "|   Moller-Plesset Perturbation Theory                                         |\n"
    @output "|       module written by M.M. Davis                                           |\n"
    @output "================================================================================\n" 
end

include("RMP2.jl")
include("UMP2.jl")
#include("DirectRMP2.jl")
include("DF-RMP2.jl")

end #module
