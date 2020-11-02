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
        throw(Fermi.InvalidFermiOption("Invalid precision: $(A)"))
    end
end

function select_algorithm(A::String)
    implemented = Dict{Any,Any}("sparse" => SparseHamiltonian(),
                                "aci" => ACI(),
                                "xaci" => xACI(),
                                "olsen" => OlsenMethod())
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOption("Invalid CI algorithm: $(A)"))
    end
end

# Implementations
include("DetOperations.jl")
include("MatrixElement.jl")
include("SparseHamiltonian.jl")
include("StringHamiltonian.jl")
include("OlsenMethod.jl")
include("ACI.jl")
include("xACI.jl")

# Most general function. It defines the precision and call a precision-specific function
"""
    Fermi.CoupledCluster.RCCSD()

Compute a RCCSD wave function using Fermi.CurrentOptions data.
"""
function CASCI(args...)
    @output "Selecting CAS precision...\n"
    prec = select_precision(Fermi.CurrentOptions["precision"])
    alg = select_algorithm(Fermi.CurrentOptions["ci_alg"])
    CASCI{prec}(args...,alg)
end

function CASCI{T}() where T <: AbstractFloat
    alg = select_algorithm(Fermi.CurrentOptions["ci_alg"])
    CASCI{T}(alg)
end

#function CASCI{T}(args...) where T <: AbstractFloat
#    @output "Selecting CAS algorithm...\n"
#    for x in args
#        print(typeof(x))
#        print("  ")
#    end
#    println()
#    alg = select_algorithm(Fermi.CurrentOptions["ci_alg"])
#    if any(i->typeof(i) <: CIAlgorithm, args)
#        throw(Fermi.MethodArgument("Invalid arguments passed to the CASCI method. Please notice that the order matters. Check docummentation for help. You can do it interactively using `? Fermi.ConfigurationInteraction.CASCI`"))
#    end
#    CASCI{T}(args...,alg)
#end
#function CASCI(args...;alg=select_algorithm(Fermi.CurrentOptions["ci_alg"]))
#    @output "Selecting CAS precision...\n"
#    prec = select_precision(Fermi.CurrentOptions["precision"])
#    CASCI{prec}(args...,alg)
#end
#
#function CASCI{T}(args...;alg=select_algorithm(Fermi.CurrentOptions["ci_alg"]))where T <: AbstractFloat
#    @output "i a iiiii\n"
#    CASCI{T}(args...,alg)
#end

"""
    Fermi.CoupledCluster.RCCSD{T}()

Compute a RCCSD wave function for a given precision T (Float64 or Float32)
"""
#function CASCI{T}(args...) where T <: AbstractFloat
#    if length(args) > 0 
#        @assert typeof(args[end]) != CIAlgorithm
#    end
#    @output "Selecting CAS Algorithm...\n"
#    alg = select_algorithm(Fermi.CurrentOptions["ci_alg"])
#    CASCI{T}(args...,alg)
#end

