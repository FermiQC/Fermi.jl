"""
    Fermi Quantum Chemistry Module

Authors: G. J. R. Aroeira, M. M. Davis, J. M. Turney, and H. F. Schaefer

GitHub: [Fermi.jl](https://github.com/FermiQC/Fermi.jl)
"""
module Fermi

include("Backend/Error.jl")
include("Core/Options.jl")                             
include("Backend/Arrays.jl")
include("Backend/Contract.jl")
include("Backend/Output.jl")
include("Core/DIIS.jl")
include("Core/PhysicalConstants.jl")
include("Core/Geometry.jl")                               
include("Core/Orbitals.jl")
include("Core/AuxiliarStructs.jl")                  
include("Backend/IntegralHelper.jl")
include("Methods/HartreeFock/HartreeFock.jl")

include("Core/SinglePointEnergy.jl")

end # module
