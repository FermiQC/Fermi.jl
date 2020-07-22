using TensorOperations

# Define Algorithims
abstract type RHFAlgorithm end
struct ConventionalRHF <: RHFAlgorithm end

"""
    Fermi.HartreeFock.RHF

Wave function object for Restricted Hartree-Fock methods

# Fields:

    molecule    Molecule object used to compute the RHF wave function
    energy      RHF Energy
    ndocc       Number of doubly occupied spatial orbitals
    nvir        Number of virtual spatial orbitals
    C           Array with MO coefficients
    eps         Array with MO energies

_struct tree:_

**RHF** <: AbstractHFWavefunction <: AbstractReferenceWavefunction <: AbstractWavefunction
"""
struct RHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    C::Array{Float64,2} 
    eps::Array{Float64,1}
end

# Algorithm-specific dispatches
include("ConventionalRHF.jl")

function select_algorithm(A::String)
    implemented = Dict{String,Any}(
        "conventional" => ConventionalRHF()
       )

    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(A)"))
    end
end

# General dispatches
"""
    Fermi.HartreeFock.RHF()

Compute RHF wave function using data from Fermi.CurrentOptions
"""
function RHF()
    molecule = Molecule()
    aoint = ConventionalAOIntegrals(molecule) 
    Alg = select_algorithm(Fermi.CurrentOptions["scf_alg"])
    RHF(molecule, aoint, Alg)
end
"""
Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals)

Compute RHF wave function using the given Molecule and Integral objects.
"""
function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals)
    Alg = select_algorithm(Fermi.CurrentOptions["scf_alg"])
    RHF(molecule, aoint, Alg)
end

"""
    Fermi.HartreeFock.RHF(wfn::RHF)

    Compute RHF wave function using the input RHF wave function (wfn) to generate a guess for orbitals.
"""
function RHF(wfn::RHF)

    Alg = select_algorithm(Fermi.CurrentOptions["scf_alg"])
    RHF(wfn, Alg)
end
