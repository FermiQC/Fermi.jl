module Fermi
import Lints
import DistributedArrays

include("Backend/Tensors.jl")                             # Top level scope

include("Backend/Environment.jl")
include("Backend/IO/Output.jl")
include("Backend/contract.jl")                            # Top level scope
include("Core/ReferenceWavefunction.jl")                  # Top level scope
include("Core/Atom.jl")                                   # Top level scope
include("Core/Molecule.jl")                               # Top level scope
include("Core/Basis.jl")                                  # Top level scope
include("Methods/HartreeFock/HartreeFock.jl")
#include("Backend/IO/Input.jl")
include("Integrals/IntegralTransformation.jl")
#include("Methods/ConfigurationInteraction/ConfigurationInteraction.jl")
#include("Methods/MollerPlesset/MollerPlesset.jl")
#include("Methods/CoupledCluster/CoupledCluster.jl")

end # module
