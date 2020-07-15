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
    nocc::Int
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
        return implemented["conventional"]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(A)"))
    end
end

# General dispatches
"""
    Fermi.HartreeFock.RHF(molecule::Molecule = Molecule(), aoint::AbstractAOIntegrals = ConventionalAOIntegrals(), Alg::RHFAlgorithm = select_algorithm(Fermi.CurrentOptions["scf_algorithm"]))

Compute RHF wave function given a Molecule, Integrals and Algorithm objects. By default data on Fermi.CurrentOptions is used.
"""
function RHF(molecule::Molecule = Molecule(), aoint::AbstractAOIntegrals = ConventionalAOIntegrals(), Alg::RHFAlgorithm = select_algorithm(Fermi.CurrentOptions["scf_algorithm"]))
    RHF(molecule, aoint, Alg)
end

"""
    Fermi.HartreeFock.RHF(wfn::RHF, Alg::RHFAlgorithm = select_algorithm(Fermi.CurrentOptions["scf_algorithm"]))

Compute RHF wave function given a RHF and Algorithim objects. The RHF wavefunction will be used to generate an initial guess of orbitals. Integrals are generated using the molecule from wfn and basis set from Fermi.CurrentOptions
"""
function RHF(wfn::RHF, Alg::RHFAlgorithm = select_algorithm(Fermi.CurrentOptions["scf_algorithm"]))
    RHF(wfn, Alg)
end
