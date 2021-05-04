"""
    Fermi Quantum Chemistry Module

Authors: G. J. R. Aroeira, M. M. Davis, J. M. Turney, and H. F. Schaefer

GitHub: [Fermi.jl](https://github.com/FermiQC/Fermi.jl)
"""
module Fermi

"""
    Fermi.AbstractWavefunction

Abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

include("Backend/Error.jl")
include("Backend/Options.jl")                             
include("Backend/Arrays.jl")
include("Backend/Contract.jl")
include("Backend/Output.jl")
include("Backend/Libcint.jl")
include("Core/DIIS.jl")
include("Core/PhysicalConstants.jl")
include("Core/Geometry.jl")                               
include("Core/BasisSet.jl")
include("Core/Orbitals.jl")
include("Core/IntegralHelper.jl")
include("Methods/HartreeFock/HartreeFock.jl")
include("Methods/MollerPlesset/MollerPlesset.jl")
include("Methods/CoupledCluster/CoupledCluster.jl")

include("Tools/SinglePointEnergy.jl")

end # module
