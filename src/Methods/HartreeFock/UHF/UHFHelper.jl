
function UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc)
    # Calculate energy
    @tensoropt A[:] := 0.5 * (H[i,j] * (Dα[i,j] + Dβ[i,j]) + Fα[i,j] * Dα[i,j] + Fβ[i,j] * Dβ[i,j])
    Ee = A[1]
    return(Ee+Vnuc)
end

function buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    # Calculate Fock matrix
    ERI = ints["ERI"]
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

function odadamping(D, Ds, F, Fs)
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
    return damp
end