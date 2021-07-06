

function UHF_core_guess(ints::IntegralHelper)
    output("Using Core Guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    F̃ = Λ * F * Λ'
    _, C̃ = diagonalize(F̃, hermitian=true)
    C = Λ * C̃
    return C, Λ
end

function UHF_gwh_guess(ints::IntegralHelper)
    # Form GWH guess
    output("Using GWH Guess")
    molecule = ints.molecule
    S = ints["S"]
    d, U = diagonalize(S, sortby = x->1/abs(x))
    Λ = FermiMDArray(Hermitian(S.data)^(-1/2))
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    H = real.(ints["T"] + ints["V"])
    nocc = molecule.Nα
    nvir = size(S,1) - nocc
    F = similar(S)

    for i = 1:nocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:nocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ         
end

function UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc, m)
    # Calculate energy
    Ee = 0
    for i in 1:m
        for j in 1:m
            Ee += 0.5 * (H[i,j]*(Dα[j,i]+Dβ[j,i]) + Fα[i,j]*Dα[j,i] + Fβ[i,j]*Dβ[j,i])
        end
    end
    return (Ee + Vnuc)
end

function buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ERI)
    # Calculate Fock matrix
    Fα .= H
    Fβ .= H
    calcJ!(Jα, Dα, ERI)
    calcJ!(Jβ, Dβ, ERI)
    calcK!(Kα, Dα, ERI)
    calcK!(Kβ, Dβ, ERI)
    Fα .+= Jα - Kα + Jβ
    Fβ .+= Jβ - Kβ + Jα
end

function calcJ!(J, D, ERI)
    # Calculate Coloumb integrals contracted with D
    @tensoropt J[i,j] = ERI[i,j,k,l] * D[l,k]
end

function calcK!(K, D, ERI)
    # Calculate Exchange integrals contracted with D
    @tensoropt K[i,j] = ERI[i,l,k,j] * D[l,k]
end

function buildD!(D, C, N)
    # Build density matrix
    Co = C[:,1:N]
    @tensoropt D[μ, ν] = Co[μ, i] * Co[ν, i]
end

function odadamping!(diis, damp, D, Ds, F, Fs)
    diis = false
    dD = D - Ds
    s = tr(Fs * dD)
    c = tr((F - Fs) * (dD))
    if c <= -s/(2*c)
        λ = 1.0
    else
        λ = -s/(2*c)
    end
    Fs .= (1-λ)*Fs + λ*F
    Ds .= (1-λ)*Ds + λ*D
    damp = 1-λ
    F .= Fs
end