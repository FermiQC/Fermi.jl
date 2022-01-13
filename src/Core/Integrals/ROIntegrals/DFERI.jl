function compute_ERI!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bμν = aoints["ERI"]

    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt (Q=>50x, p=>x, q=>x, μ=>x, ν=>x) begin
        Bpq[Q,p,q] :=  Bμν[Q, μ, ν]*C[μ, p]*C[ν, q]
    end
    I["ERI"] = Bpq
end

function compute_BOO!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bμν = aoints["ERI"]

    core = Options.get("drop_occ")
    o = (1+core):I.molecule.Nα
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Co = T.(Co)
    end
    @tensoropt (Q=>50x, i=>x, j=>x, μ=>10x, ν=>10x) begin
        Bij[Q,i,j] :=  Bμν[Q, μ, ν]*Co[μ, i]*Co[ν, j]
    end
    I["BOO"] = Bij
end

function compute_BOV!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bμν = aoints["ERI"]

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
    @tensoropt (Q=>50x, a=>8x, i=>x, μ=>10x, ν=>10x) begin
        Bia[Q,i,a] :=  Bμν[Q, μ, ν]*Co[μ, i]*Cv[ν, a]
    end
    I["BOV"] = Bia
end

function compute_BVV!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bμν = aoints["ERI"]

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
    end
    @tensoropt (Q=>50x, a=>8x, b=>8x, μ=>10x, ν=>10x) begin
        Bab[Q,a,b] :=  Bμν[Q, μ, ν]*Cv[μ, a]*Cv[ν, b]
    end
    I["BVV"] = Bab
end

function compute_OOOO!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BOO
    if haskey(I.cache, "BOO")
        Boo = I["BOO"]
    else
        compute_BOO!(I, aoints)
        Boo = I["BOO"]
        delete!(I, "BOO")
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOOO[i,j,k,l] :=  Boo[Q, i, j]*Boo[Q, k, l]
    end
    I["OOOO"] = OOOO
end

function compute_OOOV!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BOO
    if haskey(I.cache, "BOO")
        Boo = I["BOO"]
    else
        compute_BOO!(I, aoints)
        Boo = I["BOO"]
        delete!(I, "BOO")
    end

    # Get BOV
    if haskey(I.cache, "BOV")
        Bov = I["BOV"]
    else
        compute_BOV!(I, aoints)
        Bov = I["BOV"]
        delete!(I, "BOV")
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOOV[i,j,k,a] :=  Boo[Q, i, j]*Bov[Q, k, a]
    end
    I["OOOV"] = OOOV
end

function compute_OOVV!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BOO
    if haskey(I.cache, "BOO")
        Boo = I["BOO"]
    else
        compute_BOO!(I, aoints)
        Boo = I["BOO"]
        delete!(I, "BOO")
    end

    # Get BVV
    if haskey(I.cache, "BVV")
        Bvv = I["BVV"]
    else
        compute_BVV!(I, aoints)
        Bvv = I["BVV"]
        delete!(I, "BVV")
    end
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOVV[i,j,a,b] :=  Boo[Q, i, j]*Bvv[Q, a, b]
    end
    I["OOVV"] = OOVV
end

function compute_OVOV!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BOV
    if haskey(I.cache, "BOV")
        Bov = I["BOV"]
    else
        compute_BOV!(I, aoints)
        Bov = I["BOV"]
        delete!(I, "BOV")
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OVOV[i,a,j,b] :=  Bov[Q, i, a]*Bov[Q, j, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BOV
    if haskey(I.cache, "BOV")
        Bov = I["BOV"]
    else
        compute_BOV!(I, aoints)
        Bov = I["BOV"]
        delete!(I, "BOV")
    end
    
    # Get BVV
    if haskey(I.cache, "BVV")
        Bvv = I["BVV"]
    else
        compute_BVV!(I, aoints)
        Bvv = I["BVV"]
        delete!(I, "BVV")
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OVVV[i,a,b,c] :=  Bov[Q, i, a]*Bvv[Q, b, c]
    end
    I["OVVV"] = OVVV
end

function compute_VVVV!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    # Get BVV
    if haskey(I.cache, "BVV")
        Bvv = I["BVV"]
    else
        compute_BVV!(I, aoints)
        Bvv = I["BVV"]
        delete!(I, "BVV")
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        VVVV[a,b,c,d] :=  Bvv[Q, a, b]*Bvv[Q, c, d]
    end
    I["VVVV"] = VVVV
end