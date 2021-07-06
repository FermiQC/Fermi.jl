

function UHF_core_guess(ints::IntegralHelper)
    output("Using core guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    F̃ = Λ * F * Λ'
    _, C̃ = diagonalize(F̃, hermitian=true)
    C = Λ * C̃
    return C, Λ
end

function UHF_gwh_guess(ints::IntegralHelper)
    return 0
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