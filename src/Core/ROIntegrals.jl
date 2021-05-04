function mo_from_ao!(I::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}, entries...) where {T1<:AbstractFloat,T2<:AbstractFloat,
                                                            E1<:AbstractERI,E2<:AbstractERI,O<:AbstractRestrictedOrbitals}
    if T1 !== T2 || E1 !== E2 
        output("!! AO Integrals are not the same type as the MO. New integrals will be computed.")
        for entry in entries
            if occursin(r"F[dijab]{0,2}", entry)
                output("Fock matrix will be computed using old ERI")
                compute_F(I, aoints)
            end
            entries = [entries...]
            filter!(i->i==entry, entries)
        end
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    end
    t = @elapsed begin
        output("Computing MO Integrals...")
        for entry in entries
            compute!(I, entry, aoints)
        end
    end
    output("Done in {:5.5f} seconds.", t)
end

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
    elseif occursin(r"F[dijab]{0,2}", entry)
        compute_F(I, x...)
    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute!(I::IntegralHelper{T,E,O}, entry::String, x...) where {T<: AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
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
    elseif occursin(r"F[dijab]{0,2}", entry)
        compute_F(I, x...)
    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_S!(I, aoints)
end

function compute_T!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_T!(I,aoints)
end

function compute_V!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_V!(I,aoints)
end

function compute_ERI!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_ERI!(I,aoints)
end

function compute_OOOO!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Boo = I["BOO"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOOO[i,j,k,l] :=  Boo[Q, i, j]*Boo[Q, k, l]
    end
    I["OOOO"] = OOOO
end

function compute_OOOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Boo = I["BOO"]
    Bov = I["BOV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOOV[i,j,k,a] :=  Boo[Q, i, j]*Bov[Q, k, a]
    end
    I["OOOV"] = OOOV
end

function compute_OOVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Boo = I["BOO"]
    Bvv = I["BVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OOVV[i,j,a,b] :=  Boo[Q, i, j]*Bvv[Q, a, b]
    end
    I["OOVV"] = OOVV
end

function compute_OVOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Bov = I["BOV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OVOV[i,a,j,b] :=  Bov[Q, i, a]*Bov[Q, j, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Bov = I["BOV"]
    Bvv = I["BVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        OVVV[i,a,b,c] :=  Bov[Q, i, a]*Bvv[Q, b, c]
    end
    I["OVVV"] = OVVV
end

function compute_VVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    Bvv = I["BVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>20x) begin 
        VVVV[a,b,c,d] :=  Bvv[Q, a, b]*Bvv[Q, c, d]
    end
    I["VVVV"] = VVVV
end

function compute_OOOO!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_OOOO!(I,aoints)
end

function compute_OOOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_OOOV!(I,aoints)
end

function compute_OOVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_OOVV!(I,aoints)
end

function compute_OVOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_OVOV!(I,aoints)
end

function compute_OVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        compute_OVVV!(I,aoints)
end
function compute_VVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
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
        @tensoropt (Q=>50x, p=>x, q=>x, μ=>x, ν=>x) begin
            Bpq[Q,p,q] :=  Bμν[Q, μ, ν]*C[μ, p]*C[ν, q]
        end
        I["ERI"] = Bpq
end

function compute_BOO!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_BOO!(I,aoints)
end

function compute_BVV!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_BVV!(I,aoints)
end

function compute_BOV!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_BOV!(I,aoints)
end

function compute_BOO!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        core = Options.get("drop_occ")
        o = (1+core):I.molecule.Nα
        Co = I.orbitals.C[:,o]
        @tensoropt (Q=>50x, i=>x, j=>x, μ=>10x, ν=>10x) begin
            Bij[Q,i,j] :=  Bμν[Q, μ, ν]*Co[μ, i]*Co[ν, j]
        end
        I["BOO"] = Bij
end

function compute_BVV!(I::IntegralHelper{T, E, O}, aoints::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

        Bμν = aoints["ERI"]

        inac = Options.get("drop_vir")
        ndocc = I.molecule.Nα
        nbf = size(I.orbitals.C,1)
        v = (ndocc+1):(nbf - inac)
        Cv = I.orbitals.C[:,v]
        @tensoropt (Q=>50x, a=>8x, b=>8x, μ=>10x, ν=>10x) begin
            Bab[Q,a,b] :=  Bμν[Q, μ, ν]*Cv[μ, a]*Cv[ν, b]
        end
        I["BVV"] = Bab
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
        @tensoropt (Q=>50x, a=>8x, i=>x, μ=>10x, ν=>10x) begin
            Bia[Q,i,a] :=  Bμν[Q, μ, ν]*Co[μ, i]*Cv[ν, a]
        end
        I["BOV"] = Bia
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
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
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        VVVV[a,b,c,d] :=  AOERI[μ, ν, ρ, σ]*Cv[μ, a]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end

function compute_F(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_F(I, aoints)
end

function compute_F(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T, Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    ERI = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = FermiMDArray(diag(Fmol))
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F(I::IntegralHelper{T,E1,O}, aoints::IntegralHelper{T, E2, AtomicOrbitals}) where {T<:AbstractFloat, 
                                                    E1<:AbstractDFERI, E2<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    b = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = FermiMDArray(diag(Fmol))
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F(I::IntegralHelper{T,E,RHFOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    No = length(o)
    Nv = length(v)

    I["Fd"] = I.orbitals.eps
    I["Fii"] = I.orbitals.eps[o]
    I["Faa"] = I.orbitals.eps[v]
    I["Fia"] = FermiMDzeros(T, No,Nv)
    I["Fij"] = FermiMDzeros(T, No,No)
    I["Fab"] = FermiMDzeros(T, Nv,Nv)
end

function compute_ref_energy(I::IntegralHelper{T,Chonky,O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    ndocc = I.molecule.Nα
    o = 1:ndocc

    H = I["T"] + I["V"]

    OOOO = I["OOOO"]

    E0 = zero(T)
    for i in o
        E0 += 2*H[i,i]
        for j in o
            E0 += 2*OOOO[i,i,j,j] - OOOO[i,j,i,j]
        end
    end

    return E0
end

function compute_ref_energy(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat,E<:AbstractDFERI,O<:AbstractRestrictedOrbitals}
    ndocc = I.molecule.Nα
    o = 1:ndocc

    H = I["T"] + I["V"]

    Boo = I["BOO"]

    E0 = zero(T)
    for i in o
        E0 += 2*H[i,i]
        for j in o
            E0 += 2*sum(Boo[:,i,i] .* Boo[:,j,j])
            E0 -= sum(Boo[:,i,j] .* Boo[:,i,j])
        end
    end

    return E0
end