function compute!(I::IntegralHelper, entry::String)

    if entry == "S" 
        compute_S!(I)
    elseif entry == "T" 
        compute_T!(I)
    elseif entry == "V" 
        compute_V!(I)
    elseif entry == "ERI" 
        compute_ERI!(I)
    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
        I.cache["S"] = FermiMDArray(ao_overlap(I.molecule, I.basis, normalize = I.normalize))
end

function compute_T!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
        I.cache["T"] = FermiMDArray(ao_kinetic(I.molecule, I.basis, normalize = I.normalize))
end

function compute_V!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
        I.cache["V"] = FermiMDArray(ao_nuclear(I.molecule, I.basis, normalize = I.normalize))
end

function compute_ERI!(I::IntegralHelper{T, JKFIT, AtomicOrbitals}) where T<:AbstractFloat
        auxjk = Options.get("jkfit")
        basis = I.basis

        # If aux is auto, determine the aux basis from the basis
        if auxjk == "auto"
            std_name = Regex("cc-pv.z")
            auxjk = occursin(std_name, basis) ? basis*"-jkfit" : "cc-pvqz-jkfit"
        end

        I.cache["ERI"] = FermiMDArray(df_ao_eri(I.molecule, I.basis, auxjk, normalize = I.normalize))
end

function compute_ERI!(I::IntegralHelper{T, RIFIT, AtomicOrbitals}) where T<:AbstractFloat
        auxri = Options.get("rifit")
        basis = I.basis

        # If aux is auto, determine the aux basis from the basis
        if auxri == "auto"
            std_name = Regex("cc-pv.z")
            auxri = occursin(std_name, basis) ? basis*"-rifit" : "cc-pvqz-rifit"
        end

        I.cache["ERI"] = FermiMDArray(df_ao_eri(I.molecule, I.basis, auxri, normalize = I.normalize))
end

function compute_ERI!(I::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat
        I.cache["ERI"] = FermiMDArray(ao_eri(I.molecule, I.basis, normalize = I.normalize))
end