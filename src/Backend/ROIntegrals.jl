function compute!(I::IntegralHelper{T,Chonky,O}, entry::String, x...) where {T<: AbstractFloat,O<:AbstractRestrictedOrbitals}
    if entry == "S"
        compute_S!(I, x...)
    elseif entry == "T"
        compute_T!(I, x...)
    elseif entry == "V"
        compute_V!(I, x...)
    elseif entry == "ERI"
        compute_ERI!(I, x...)
    elseif entry == "OOOO"
        compute_OOOO!(I, x...)
    elseif entry == "OOOV"
        compute_OOOV!(I, x...)
    elseif entry == "OVOV"
        compute_OVOV!(I, x...)
    elseif entry == "OOVV"
        compute_OOVV!(I, x...)
    elseif entry == "OVVV"
        compute_OVVV!(I, x...)
    elseif entry == "VVVV"
        compute_VVVV!(I, x...)
    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute!(I::IntegralHelper{T,RIFIT,O}, entry::String, x...) where {T<: AbstractFloat,O<:AbstractRestrictedOrbitals}
    if entry == "S"
        compute_S!(I, x...)
    elseif entry == "T"
        compute_T!(I, x...)
    elseif entry == "V"
        compute_V!(I, x...)
    elseif entry == "ERI"
        compute_ERI!(I, x...)
    elseif entry == "BOO"
        compute_BOO!(I, x...)
    elseif entry == "BOV"
        compute_BOV!(I, x...)
    elseif entry == "BVV"
        compute_BVV!(I, x...)
    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_S!(I, aoints)
end

function compute_T!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_T!(I,aoints)
end

function compute_V!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_V!(I,aoints)
end

function compute_ERI!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_ERI!(I,aoints)
end

function compute_OOOO!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_OOOO!(I,aoints)
end

function compute_OOOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_OOOV!(I,aoints)
end

function compute_OOVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_OOVV!(I,aoints)
end

function compute_OVOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_OVOV!(I,aoints)
end

function compute_OVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_OVVV!(I,aoints)
end
function compute_VVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_VVVV!(I,aoints)
end

function compute_S!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Compute AO overlap
        Sao = ints["S"]

        # Convert AO overlap to MO overlap
        C = I.orbitals.C
        @tensoropt S[i,j] := Sao[μ, ν]*C[μ,i]*C[ν,j]
        I.cache["S"] = S
end

function compute_T!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Compute AO overlap
        Tao = ints["T"]

        # Convert AO overlap to MO overlap
        C = I.orbitals.C
        @tensoropt Tmo[i,j] := Tao[μ, ν]*C[μ,i]*C[ν,j]
        I.cache["T"] = Tmo
end

function compute_V!(I::IntegralHelper{T, E, O}, ints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Compute AO overlap
        Vao = ints["V"]

        # Convert AO overlap to MO overlap
        C = I.orbitals.C
        @tensoropt V[i,j] := Vao[μ, ν]*C[μ,i]*C[ν,j]
        I.cache["V"] = V
end

function compute_ERI!(I::IntegralHelper{T, Chonky, O}, aoints::IntegralHelper{T, Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
        AOERI = aoints["ERI"]
        C = I.orbitals.C
        @tensoropt MOERI[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*C[μ, i]*C[ν, j]*C[ρ, k]*C[σ, l]
        I["ERI"] = MOERI
end

function compute_ERI!(I::IntegralHelper{T, RIFIT, O}, aoints::IntegralHelper{T, RIFIT, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        C = I.orbitals.C
        @tensoropt (Q=>5, p=>1, q=>1, μ=>1, ν=>1) begin
            Bpq[Q,p,q] :=  Bμν[Q, μ, ν]*C[μ, p]*C[ν, q]
        end
        I["ERI"] = Bpq
end

function compute_BOO!(I::IntegralHelper{T, RIFIT, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_BOO!(I,aoints)
end

function compute_BVV!(I::IntegralHelper{T, RIFIT, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_BVV!(I,aoints)
end

function compute_BOV!(I::IntegralHelper{T, RIFIT, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        aoints = IntegralHelper(I, AtomicOrbitals())
        compute_BOV!(I,aoints)
end

function compute_BOO!(I::IntegralHelper{T, RIFIT, O}, aoints::IntegralHelper{T, RIFIT, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        core = Options.get("drop_occ")
        o = (1+core):I.molecule.Nα
        Co = I.orbitals.C[:,o]
        @tensoropt (Q=>20, i=>1, j=>1, μ=>10, ν=>10) begin
            Bij[Q,i,j] :=  Bμν[Q, μ, ν]*Co[μ, i]*Co[ν, j]
        end
        I["BOO"] = Bij
end

function compute_BVV!(I::IntegralHelper{T, RIFIT, O}, aoints::IntegralHelper{T, RIFIT, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        inac = Options.get("drop_vir")
        ndocc = I.molecule.Nα
        nvir = size(I.orbitals.C,1) - ndocc
        v = (ndocc+1):(nvir - inac)
        Cv = I.orbitals.C[:,v]
        @tensoropt (Q=>50, a=>8, b=>8, μ=>10, ν=>10) begin
            Bab[Q,a,b] :=  Bμν[Q, μ, ν]*Cv[μ, a]*Cv[ν, b]
        end
        I["BVV"] = Bab
end

function compute_BOV!(I::IntegralHelper{T, RIFIT, O}, aoints::IntegralHelper{T, RIFIT, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        inac = Options.get("drop_vir")
        core = Options.get("drop_occ")
        ndocc = I.molecule.Nα
        nvir = size(I.orbitals.C,1) - ndocc
        o = (1+core):ndocc
        v = (ndocc+1):(nvir - inac)
        Cv = I.orbitals.C[:,v]
        Co = I.orbitals.C[:,o]
        @tensoropt (Q=>50, a=>8, i=>1, μ=>10, ν=>10) begin
            Bia[Q,i,a] :=  Bμν[Q, μ, ν]*Co[μ, i]*Cv[ν, a]
        end
        I["BOV"] = Bia
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        OOOO[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
    end
    I["OOOO"] = OOOO
end

function compute_OOOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nvir = size(I.orbitals.C,1) - ndocc
    o = (1+core):ndocc
    v = (ndocc+1):(nvir - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        OOOV[i,j,k,a] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
    end
    I["OOOV"] = OOOV
end

function compute_OOVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nvir = size(I.orbitals.C,1) - ndocc
    o = (1+core):ndocc
    v = (ndocc+1):(nvir - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        OOVV[i,j,a,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
    end
    I["OOVV"] = OOVV
end

function compute_OVOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nvir = size(I.orbitals.C,1) - ndocc
    o = (1+core):ndocc
    v = (ndocc+1):(nvir - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        OVOV[i,a,j,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nvir = size(I.orbitals.C,1) - ndocc
    o = (1+core):ndocc
    v = (ndocc+1):(nvir - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        OVVV[i,a,b,c] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
    end
    I["OVVV"] = OVVV
end

function compute_VVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nvir = size(I.orbitals.C,1) - ndocc
    v = (ndocc+1):(nvir - inac)
    Cv = I.orbitals.C[:,v]
    @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
        VVVV[a,b,c,d] :=  AOERI[μ, ν, ρ, σ]*Cv[μ, a]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end