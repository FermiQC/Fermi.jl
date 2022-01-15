function compute!(I::IntegralHelper{T,<:AbstractERI,<:AbstractRestrictedOrbitals}, entry::String) where T <: AbstractFloat
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)

    compute!(I, aoints, entry)
end

function compute!(I::IntegralHelper{T,<:AbstractERI,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,<:AbstractERI,AtomicOrbitals}, entry::String) where T<: AbstractFloat
    if entry == "S"
        compute_S!(I, aoints)
    elseif entry == "T"
        compute_T!(I, aoints)
    elseif entry == "V"
        compute_V!(I, aoints)
    elseif entry == "ERI"
        compute_ERI!(I, aoints)
    elseif entry == "OOOO"
        compute_OOOO!(I, aoints)
    elseif entry == "OOOV"
        compute_OOOV!(I, aoints)
    elseif entry == "OVOV"
        compute_OVOV!(I, aoints)
    elseif entry == "OOVV"
        compute_OOVV!(I, aoints)
    elseif entry == "OVVV"
        compute_OVVV!(I, aoints)
    elseif entry == "VVVV"
        compute_VVVV!(I, aoints)
    elseif occursin(r"F[dijab]{0,2}", entry)
        compute_F!(I, aoints)
    elseif I.eri_type === Chonky() # Out of options for Chonky
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    elseif entry == "BOO"
        compute_BOO!(I, aoints)
    elseif entry == "BOV"
        compute_BOV!(I, aoints)
    elseif entry == "BVV"
        compute_BVV!(I, aoints)
    end
end

include("OneElectron.jl")
include("DFERI.jl")
include("Chonky.jl")
include("Sparse.jl")
include("Fock.jl")
