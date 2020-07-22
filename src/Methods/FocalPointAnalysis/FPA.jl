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

    extrapolated, ext_ind = extrapolate(data, basis, methods)
    push!(basis, "CBS")

    @output "\n ⇒ Final FPA (absolute):\n"
    fpa = FPA{T}(basis, methods, extrapolated, ext_ind)
    print_fpa(fpa)
    @output "⇒ Final FPA Energy:   {:15.10f}\n" fpa.data[end,end]
    return fpa
end

function get_fpa_row(molecule::Molecule, upmethod::String, basis::String, T::DataType)

    Fermi.CurrentOptions["basis"] = basis
    
    @output "\n   → Starting FPA row for {}\n\n" basis
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
