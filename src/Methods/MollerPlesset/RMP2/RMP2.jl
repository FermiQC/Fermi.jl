using Fermi.HartreeFock: RHF

export RMP2

abstract type RMP2Algorithm end

function get_mp2_alg()
    implemented = [RMP2a()]
    N = Options.get("mp2_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(InvalidFermiOption("implementation number $N not available for RMP2."))
    end
end

"""
    RMP2{T} <: AbstractMPWavefunction

    TODO
"""
struct RMP2{T} <: AbstractMPWavefunction
    reference::RHF
    correlation::T
    energy::T
end

function RMP2(x...)
    precision = Options.get("precision")
    if precision == "single"
        RMP2{Float32}(x...)
    elseif precision == "double"
        RMP2{Float64}(x...)
    else
        throw(InvalidFermiOption("precision can only be `single` or `double`. Got $precision"))
    end
end

function RMP2{Float64}(mol::Molecule = Molecule(), ints::IntegralHelper{Float64} = IntegralHelper{Float64}())

    # Compute Restricted Hartree-Fock 
    refwfn = Fermi.HartreeFock.RHF(mol, ints)

    output("MP2 will be computed using double precision")

    # Get rid of integrals that will be no longer used
    delete!(ints, "JKERI", "S", "T", "V")

    RMP2{Float64}(refwfn, ints)
end

function RMP2{Float32}(mol::Molecule = Molecule(), ints::IntegralHelper{Float32} = IntegralHelper{Float32}())

    # Compute Restricted Hartree-Fock 
    refwfn = Fermi.HartreeFock.RHF(mol)

    output("MP2 will be computed using single precision")

    RMP2{Float32}(refwfn, ints)
end

function RMP2{T}(refwfn::RHF, ints::IntegralHelper{T} = IntegralHelper{T}()) where T <: AbstractFloat

    # MP2 only needs a reference and ERI
    if Options.get("df") 
        output("Computing RI density fitted ERI")
        t = @elapsed ints["RIERI"]
        output("Done in {:5.5f} s", t)
        RMP2{T}(refwfn, ints, get_mp2_alg())
    else
        output("Computing ERI")
        t = @elapsed ints["ERI"]
        output("Done in {:5.5f} s", t)
        RMP2{T}(refwfn, ints, get_mp2_alg())
    end
end

# For each implementation a singleton type must be create
struct RMP2a <: RMP2Algorithm end
include("RMP2a.jl")