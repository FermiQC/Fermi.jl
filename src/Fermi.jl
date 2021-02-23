"""
    Fermi Quantum Chemistry Module

Authors: M. M. Davis, G. J. R. Aroeira, J. M. Turney, and H. F. Schaefer

GitHub: [Fermi.jl](https://github.com/FermiQC/Fermi.jl)
"""
module Fermi
import Lints
import TBLIS
import DistributedArrays

include("Core/Options.jl")                             # Top level scope
include("Backend/Tensors.jl")                             # Top level scope
include("Backend/Environment.jl")
include("Backend/ComputeEnvironment.jl")
include("Backend/IO/Output.jl")
include("Backend/contract.jl")                            
include("Backend/DIIS.jl")
include("Core/PhysicalConstants.jl")
include("Core/AbstractWavefunctions.jl")                  
include("Core/Geometry.jl")                               
include("Core/Orbitals.jl")
include("Core/Integrals.jl")
include("Methods/HartreeFock/HartreeFock.jl")
include("Methods/MollerPlesset/MollerPlesset.jl")
include("Methods/ConfigurationInteraction/ConfigurationInteraction.jl")
include("Methods/CoupledCluster/CoupledCluster.jl")
#include("Methods/FocalPointAnalysis/FocalPointAnalysis.jl")

include("Core/SinglePointEnergy.jl")

end # module
