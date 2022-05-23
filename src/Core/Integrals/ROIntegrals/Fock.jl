function compute_F!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]
    if eltype(I.orbitals.C) !== T
        Co = T.(Cv)
    end

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    ERI = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = diag(Fmol)
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F!(I::IntegralHelper{T1,<:AbstractDFERI,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T2,<:AbstractDFERI,AtomicOrbitals}) where {T1<:AbstractFloat,T2<:AbstractFloat}

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]
    if eltype(I.orbitals.C) !== T1
        Co = T1.(Co)
    end

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    b = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = diag(Fmol)
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F!(I::IntegralHelper{<:AbstractFloat,Chonky,RHFOrbitals}, aoints::IntegralHelper{<:AbstractFloat, SparseERI, AtomicOrbitals})
    compute_F!(I)
end

function compute_F!(I::IntegralHelper{<:AbstractFloat,RIFIT,RHFOrbitals}, aoints::IntegralHelper{<:AbstractFloat, <:AbstractDFERI, AtomicOrbitals})
    compute_F!(I)
end

function compute_F!(I::IntegralHelper{T, <:AbstractERI, RHFOrbitals}) where T<:AbstractFloat

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    No = length(o)
    Nv = length(v)

    I["Fd"] = T.(I.orbitals.eps)
    I["Fii"] = T.(I.orbitals.eps[o])
    I["Faa"] = T.(I.orbitals.eps[v])
    I["Fia"] = zeros(T, No,Nv)
    I["Fij"] = zeros(T, No,No)
    I["Fab"] = zeros(T, Nv,Nv)
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