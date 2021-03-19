import Fermi.HartreeFock: RHF

struct RMP2{T} <: AbstractMPWavefunction# where T <: AbstractFloat
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
        throw(InvalidFermiOption("precision can only be single or double. Got $precision"))
    end
end

function RMP2{Float64}()

    output("MP2 will be computed using double precision")

    # Compute Restricted Hartree-Fock 
    mol = Molecule()
    ints = IntegralHelper(mol)
    refwfn = Fermi.HartreeFock.RHF(mol, ints)

    # Get rid of integrals that will be no longer used
    delete!(ints, "JKERI", "S", "T", "V")

    # MP2 only needs a reference and ERI
    mp2_type = Options.get("mp2_type")
    if mp2_type == "df" 

        output("Computing RI density fitted ERI")
        t = @elapsed ints["RIERI"]
        output("Done in {:5.5f} s", t)
        RMP2(refwfn, ints)

    elseif mp2_type == "conventional"
        output("Computing ERI")
        t = @elapsed ints["ERI"]
        output("Done in {:5.5f} s", t)
        RMP2(refwfn, ints)
    end
end

function RMP2{Float32}()

    output("MP2 will be computed using single precision")

    # Compute Restricted Hartree-Fock 
    refwfn = Fermi.HartreeFock.RHF()

    # Get new integrals object
    ints = IntegralHelper{Float32}()

    # MP2 only needs a reference and ERI
    mp2_type = Options.get("mp2_type")
    if mp2_type == "df" 

        output("Computing RI density fitted ERI")
        t = @elapsed ints["RIERI"]
        output("Done in {:5.5f} s", t)
        RMP2(refwfn, ints)

    elseif mp2_type == "conventional"
        output("Computing ERI")
        t = @elapsed ints["ERI"]
        output("Done in {:5.5f} s", t)
        RMP2(refwfn, ints)
    end
end

include("ConventionalRMP2.jl")
