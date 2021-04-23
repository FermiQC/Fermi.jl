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
        if Options.get("lints")
            I.cache["S"] = FermiMDArray(ao_overlap(I.molecule, I.basis, normalize = I.normalize))
        else
            bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
            I.cache["S"] = FermiMDArray(Fermi.GaussianBasis.ao_1e(bs, "overlap"))
        end
end

function compute_T!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
        if Options.get("lints")
            I.cache["T"] = FermiMDArray(ao_kinetic(I.molecule, I.basis, normalize = I.normalize))
        else
            bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
            I.cache["T"] = FermiMDArray(Fermi.GaussianBasis.ao_1e(bs, "kinetic"))
        end
end

function compute_V!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
        if Options.get("lints")
            I.cache["V"] = FermiMDArray(ao_nuclear(I.molecule, I.basis, normalize = I.normalize))
        else
            bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
            I.cache["V"] = FermiMDArray(Fermi.GaussianBasis.ao_1e(bs, "nuclear"))
        end
end

function compute_ERI!(I::IntegralHelper{T, JKFIT, AtomicOrbitals}) where T<:AbstractFloat
    auxjk = Options.get("jkfit")
    basis = I.basis

    # If aux is auto, determine the aux basis from the basis
    if auxjk == "auto"
        std_name = Regex("cc-pv.z")
        auxjk = occursin(std_name, basis) ? basis*"-jkfit" : "cc-pvqz-jkfit"
    end

    if Options.get("lints")
        I.cache["ERI"] = FermiMDArray(df_ao_eri(I.molecule, I.basis, auxjk, normalize = I.normalize))
    else
        bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
        auxbs = Fermi.GaussianBasis.BasisSet(I.molecule, auxjk)
        J = FermiMDArray(Fermi.GaussianBasis.ao_2e2c(auxbs))
        Pqp = FermiMDArray(Fermi.GaussianBasis.ao_2e3c(bs, auxbs))
        #display(Pqp)
        Jh = Array(real(J^(-1/2)))
        @tensor b[Q,p,q] := Pqp[p,q,P]*Jh[P,Q]
        I.cache["ERI"] = b
    end
end

function compute_ERI!(I::IntegralHelper{T, RIFIT, AtomicOrbitals}) where T<:AbstractFloat
    auxri = Options.get("rifit")
    basis = I.basis

    # If aux is auto, determine the aux basis from the basis
    if auxri == "auto"
        std_name = Regex("cc-pv.z")
        auxri = occursin(std_name, basis) ? basis*"-rifit" : "cc-pvqz-rifit"
    end

    if Options.get("lints")
        I.cache["ERI"] = FermiMDArray(df_ao_eri(I.molecule, I.basis, auxri, normalize = I.normalize))
    else
        bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
        auxbs = Fermi.GaussianBasis.BasisSet(I.molecule, auxri)
        J = FermiMDArray(Fermi.GaussianBasis.ao_2e2c(auxbs))
        Pqp = FermiMDArray(Fermi.GaussianBasis.ao_2e3c(bs, auxbs))
        #display(Pqp)
        Jh = Array(real(J^(-1/2)))
        @tensor b[Q,p,q] := Pqp[p,q,P]*Jh[P,Q]
        I.cache["ERI"] = b
    end
end

function compute_ERI!(I::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat
        if Options.get("lints")
            I.cache["ERI"] = FermiMDArray(ao_eri(I.molecule, I.basis, normalize = I.normalize))
        else
            bs = Fermi.GaussianBasis.BasisSet(I.molecule, I.basis)
            I.cache["ERI"] = FermiMDArray(Fermi.GaussianBasis.ao_2e4c(bs))
        end
end