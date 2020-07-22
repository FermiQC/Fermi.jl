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
include("Core/PhysicalConstants.jl")
include("Core/AbstractWavefunctions.jl")                  
include("Core/Geometry.jl")                               
include("Core/Integrals.jl")
include("Methods/HartreeFock/HartreeFock.jl")
#include("Methods/MollerPlesset/MollerPlesset.jl")
include("Methods/ConfigurationInteraction/ConfigurationInteraction.jl")
include("Methods/CoupledCluster/CoupledCluster.jl")
include("Methods/FocalPointAnalysis/FPA.jl")

end # module
