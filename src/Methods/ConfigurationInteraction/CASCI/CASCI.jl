"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant

structure holding alpha and beta strings for a determinant as integers.

# Fields

    α integer representing an alpha string
    β integer representing an beta string
"""
struct Determinant
    α::Integer
    β::Integer
end

struct CASCI{T} <: AbstractCIWavefunction
    ref::Fermi.HartreeFock.RHF
    energy::T
    dets::Array{Determinant,1}
    coef::Array{T,1}
end

function select_precision(A::String)
    implemented = Dict{Any,Any}("single" => Float32,
                                "double" => Float64)
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid precision: $(A)"))
    end
end

function select_algorithm(A::String)
    implemented = Dict{Any,Any}("sparse" => SparseHamiltonian(),
                                "aci" => ACI())
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOptions("Invalid CI algorithm: $(A)"))
    end
end

# Implementations
include("DetOperations.jl")
include("MatrixElement.jl")
include("SparseHamiltonian.jl")
include("ACI.jl")

# Most general function. It defines the precision and call a precision-specific function
"""
    Fermi.CoupledCluster.RCCSD()

Compute a RCCSD wave function using Fermi.CurrentOptions data.
"""
function CASCI()
    @output "Selecting CAS precision...\n"
    prec = select_precision(Fermi.CurrentOptions["precision"])
    CASCI{prec}()
end

"""
    Fermi.CoupledCluster.RCCSD{T}()

Compute a RCCSD wave function for a given precision T (Float64 or Float32)
"""
function CASCI{T}() where T <: AbstractFloat
    @output "Selecting CAS Algorithm...\n"
    alg = select_algorithm(Fermi.CurrentOptions["ci_alg"])
    CASCI{T}(alg)
end

