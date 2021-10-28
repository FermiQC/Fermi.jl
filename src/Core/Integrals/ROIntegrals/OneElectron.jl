function compute_S!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_S!(I, aoints)
end

function compute_T!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_T!(I,aoints)
end

function compute_V!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_V!(I,aoints)
end

function compute_S!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Compute AO overlap
    Sao = ints["S"]

    # Convert AO overlap to MO overlap
    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt S[i,j] := Sao[μ, ν]*C[μ,i]*C[ν,j]
    I.cache["S"] = S
end

function compute_T!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Compute AO overlap
    Tao = ints["T"]

    # Convert AO overlap to MO overlap
    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt Tmo[i,j] := Tao[μ, ν]*C[μ,i]*C[ν,j]
    I.cache["T"] = Tmo
end

function compute_V!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Compute AO overlap
    Vao = ints["V"]

    # Convert AO overlap to MO overlap
    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt V[i,j] := Vao[μ, ν]*C[μ,i]*C[ν,j]
    I.cache["V"] = V
end