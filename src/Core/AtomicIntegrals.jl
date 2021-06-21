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
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["S"] = FermiMDArray(ao_1e(bs, "overlap", T))
end

function compute_T!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["T"] = FermiMDArray(ao_1e(bs, "kinetic", T))
end

function compute_V!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["V"] = FermiMDArray(ao_1e(bs, "nuclear", T))
end

function compute_ERI!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI}

    bs = I.orbitals.basisset
    auxbs = I.eri_type.basisset
    J = FermiMDArray(ao_2e2c(auxbs, T))
    Pqp = FermiMDArray(ao_2e3c(bs, auxbs, T))
    Jh = Array(real(J^(-1/2)))
    @tensor b[Q,p,q] := Pqp[p,q,P]*Jh[P,Q]
    I.cache["ERI"] = b
end

function compute_ERI!(I::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = FermiMDArray(ao_2e4c(bs, T))
end

function compute_ERI!(I::IntegralHelper{T, SparseERI, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = Fermi.FermiSparse(sparse_ao_2e4c(bs, T)...)
end

function index2(i::Signed, j::Signed)::Signed
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

# Produces all unique indices ijkl for the two-electron integral
# Note, these indexes start from zero, because we use it in a ccall
function find_indices(nbf::Signed)

    out = NTuple{4,Int16}[]
    N = Int16(nbf - 1)
    ZERO = zero(Int16)

    for i = ZERO:N
        for j = i:N
            for k = ZERO:N
                for l = k:N
                    if index2(i,j) < index2(k,l)
                        continue
                    end
                    push!(out, (i,j,k,l))
                end
            end
        end
    end

    return out
end

include("Integrals/OneElectron.jl")
include("Integrals/TwoElectronTwoCenter.jl")
include("Integrals/TwoElectronThreeCenter.jl")
include("Integrals/TwoElectronFourCenter.jl")