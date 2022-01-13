function compute!(I::IntegralHelper{<:AbstractFloat,<:AbstractERI,AtomicOrbitals}, entry::String)

    if entry == "S" 
        compute_S!(I)
    elseif entry == "T" 
        compute_T!(I)
    elseif entry == "V" 
        compute_V!(I)
    elseif entry == "ERI" 
        compute_ERI!(I)
    else
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["S"] = FermiMDArray(overlap(bs, T))
end

function compute_T!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["T"] = FermiMDArray(kinetic(bs, T))
end

function compute_V!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["V"] = FermiMDArray(nuclear(bs, T))
end

function compute_ERI!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI}

    bs = I.orbitals.basisset
    auxbs = I.eri_type.basisset
    J = FermiMDArray(ERI_2e2c(auxbs, T))
    Pqp = FermiMDArray(ERI_2e3c(bs, auxbs, T))
    Jh = similar(J)
    Jh .= Array(real(J^(-1/2)))
    @tensor b[Q,p,q] := Pqp[p,q,P]*Jh[P,Q]
    I.cache["ERI"] = b
end

function compute_ERI!(I::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = FermiMDArray(ERI_2e4c(bs, T))
end

function compute_ERI!(I::IntegralHelper{T, SparseERI, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = Fermi.FermiSparse(sparseERI_2e4c(bs, T)...)
end