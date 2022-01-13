function compute_ERI!(I::IntegralHelper{T, Chonky, O}, aoints::IntegralHelper{T, Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    AOERI = aoints["ERI"]
    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt MOERI[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*C[μ, i]*C[ν, j]*C[ρ, k]*C[σ, l]
    I["ERI"] = MOERI
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOO[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
    end
    I["OOOO"] = OOOO
end

function compute_OOOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOV[i,j,k,a] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
    end
    I["OOOV"] = OOOV
end

function compute_OOVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOVV[i,j,a,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
    end
    I["OOVV"] = OOVV
end

function compute_OVOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVOV[i,a,j,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVVV[i,a,b,c] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
    end
    I["OVVV"] = OVVV
end

function compute_VVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        VVVV[a,b,c,d] :=  AOERI[μ, ν, ρ, σ]*Cv[μ, a]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end