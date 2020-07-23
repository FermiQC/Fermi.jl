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
    @set basis 3-21g
    currentRHF = RHF(molecule)
    for (b,m) in grid
        currentRHF, row = get_fpa_row(molecule, m, b, T, currentRHF)
        push!(basis, b)
        for i in eachindex(row)
            data[x,i] = row[i]
        end
        x += 1
    end

    extrapolated, ext_ind = extrapolate(data, basis, methods)
    push!(basis, "CBS")

    fpa = FPA{T}(basis, methods, extrapolated, ext_ind)
    @output "\n⇒ Final FPA Energy:   {:15.10f}\n" fpa.data[end,end]
    print_fpa(fpa)
    return fpa
end

function get_fpa_row(molecule::Molecule, upmethod::String, basis::String, T::DataType, RHFguess::Nothing)

    Fermi.CurrentOptions["basis"] = basis
    
    @output "\n   → Starting FPA row for {}\n\n" basis
    aoint = ConventionalAOIntegrals(molecule)
    RHFwfn = RHF(molecule, aoint)

    if lowercase(upmethod) == "rhf"
        return RHFwfn, (RHFwfn.energy)
    end

    get_fpa_corr_row(upmethod, basis, T, RHFwfn, aoint)
end

function get_fpa_row(molecule::Molecule, upmethod::String, basis::String, T::DataType, RHFguess::RHF)

    Fermi.CurrentOptions["basis"] = basis
    
    @output "\n   → Starting FPA row for {}\n\n" basis
    aoint = ConventionalAOIntegrals(molecule)
    RHFwfn = RHF(RHFguess, aoint)

    if lowercase(upmethod) == "rhf"
        return RHFwfn, (RHFwfn.energy)
    end

    get_fpa_corr_row(upmethod, basis, T, RHFwfn, aoint)
end

function get_fpa_corr_row(upmethod::String, basis::String, T::DataType, RHFwfn::RHF, aoint::ConventionalAOIntegrals)

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    @output "Transforming Integrals..."
    moint = Fermi.Integrals.PhysRestrictedMOIntegrals{T}(RHFwfn.ndocc, RHFwfn.nvir, drop_occ, drop_vir, RHFwfn.C, aoint)
    aoint = nothing

    CCSDwfn = RCCSD{T}(RHFwfn, moint, Fermi.CoupledCluster.CTF()) 

    if lowercase(upmethod) == "ccsd"
        return RHFwfn, (RHFwfn.energy, CCSDwfn.GuessEnergy, CCSDwfn.CorrelationEnergy)
    end

    CCSDpTwfn = RCCSDpT{T}(CCSDwfn, moint)

    if lowercase(upmethod) == "ccsd(t)"
        return RHFwfn, (RHFwfn.energy, CCSDwfn.GuessEnergy, CCSDwfn.CorrelationEnergy, CCSDpTwfn.correction+CCSDwfn.CorrelationEnergy)
    end
end
