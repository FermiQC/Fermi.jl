module FocalPointAnalysis

using Fermi
using Fermi.Output
using Fermi.Integrals: ConventionalAOIntegrals, PhysRestrictedMOIntegrals
using Fermi.HartreeFock: RHF
using Fermi.Geometry: Molecule
using Fermi.CoupledCluster: RCCSD, RCCSDpT
using PrettyTables

struct FPA{T}
    basis::Array{String,1}
    methods::Array{String,1}
    data::Array{T,2}
    extrapolated_indexes::Array{Tuple{Int64,Int64},1}
end

function get_fpa_row(molecule::Molecule, upmethod::String, basis::String, T::DataType)

    Fermi.CurrentOptions["basis"] = basis
    
    @output "   → Starting FPA row for {}\n\n" basis
    aoint = Fermi.Integrals.ConventionalAOIntegrals(molecule)
    RHFwfn = Fermi.HartreeFock.RHF(molecule, aoint)

    if lowercase(upmethod) == "rhf"
        return RHFwfn.energy
    end

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    @output "Transforming Integrals..."
    moint = Fermi.Integrals.PhysRestrictedMOIntegrals{T}(RHFwfn.ndocc, RHFwfn.nvir, drop_occ, drop_vir, RHFwfn.C, aoint)
    aoint = nothing

    CCSDwfn = RCCSD{T}(RHFwfn, moint, Fermi.CoupledCluster.CTF()) 

    if lowercase(upmethod) == "ccsd"
        return RHFwfn.energy, CCSDwfn.GuessEnergy, CCSDwfn.CorrelationEnergy
    end

    CCSDpTwfn = RCCSDpT{T}(CCSDwfn, moint)

    if lowercase(upmethod) == "ccsd(t)"
        return RHFwfn.energy, CCSDwfn.GuessEnergy, CCSDwfn.CorrelationEnergy, CCSDpTwfn.correction+CCSDwfn.CorrelationEnergy
    end
end

function FPA(grid::Array{Tuple{String,String},1})

    if Fermi.CurrentOptions["precision"] == "single"
        prec = Float32
    elseif Fermi.CurrentOptions["precision"] == "double"
        prec = Float64
    else
        throw(Fermi.InvalidFermiOptions("Invalid precision: $(A)"))
    end

    FPA{prec}(grid)
end

function FPA{T}(grid::Array{Tuple{String,String},1}) where T <: AbstractFloat
    methods = ["RHF", "MP2", "CCSD", "CCSD(T)"]
    basis = String[] 
    data = zeros(length(grid), 4)
    molecule = Molecule()
    x = 1
    for (b,m) in grid
        row = get_fpa_row(molecule, m, b, T)
        push!(basis, b)
        for i in eachindex(row)
            data[x,i] = row[i]
        end
        x += 1
    end

    fpa = FPA{T}(basis, methods, data)

    @output "\n ⇒ Final FPA (absolute):\n"
    print_fpa(fpa)
    return fpa
end

function extrapolate(fpa::FPA)
    X = Int64[]
    for b in fpa.basis
occursin
issorted
end

function feller_formula(E::Array{T,1}, X::Array{Int64,1}) where T <: AbstractFloat
    
    # SCF extrapolation: Escf = A + B*exp(-cX)
    # From: Feller 1993 (https://doi.org/10.1063/1.464749)

    # Make sure E and X are ordered
    sp = sortperm(X)
    E = E[sp]
    X = X[sp]

    C = log((E[2] - E[1])/(E[3] - E[2]))
    B = (E[3] - E[2])/(exp(-C*X[3]) - exp(-C*X[2]))
    A = E[3] - B*exp(-C*X[3])
    
    return A,B,C
end

function helgaker_formula(E::Array{T,1}, X::Array{Int64,1}) where T <: AbstractFloat

    # Correlation energy extrapolation: Ecorr = A + B*X^-3
    # From Helgaker 1997 (https://doi.org/10.1063/1.473863)

    # Make sure E and X are ordered
    sp = sortperm(X)
    E = E[sp]
    X = X[sp]

    B = (E[2] - E[1])/(X[2]^-3 - X[1]^-3)
    A = E[2] - B*X[2]^-3
    return A,B
end

function print_fpa(fpa::FPA)
    rmv_zeros = (v,i,j) -> v != 0.0 ? v : " "

    out = pretty_table(String, fpa.data, fpa.methods, alignment=:c,
    formatters = (rmv_zeros, ft_printf("%15.10f")), row_names=fpa.basis, 
    row_name_column_title="Basis Set")

    @output "\n{}\n" out
end
end #module
