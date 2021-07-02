function UHF(Alg::UHFa)
    ints = IntegralHelper{Float64}(eri_type=Chonky())
    Fermi.HartreeFock.hf_header()
    molecule = ints.molecule
    Nα = molecule.Nα
    Nβ = molecule.Nβ
    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    G = ints["ERI"]
    X = S^(-1/2)
    m = size(S)[1]
    Dα = FermiMDzeros(Float64, (m,m))
    Dβ = FermiMDzeros(Float64, (m,m))
    Jα = FermiMDzeros(Float64, (m,m))
    Jβ = FermiMDzeros(Float64, (m,m))
    Kα = FermiMDzeros(Float64, (m,m))
    Kβ = FermiMDzeros(Float64, (m,m))
    E = 0
    Cα = FermiMDzeros(Float64, (m,m))
    Cβ = FermiMDzeros(Float64, (m,m))
    ϵα = FermiMDzeros(Float64, (m))
    ϵβ = FermiMDzeros(Float64, (m))
    conv = Options.get("scf_e_conv")
    max_iter = Options.get("scf_max_iter")
    converged = false
    ite = 0
    while ite <= max_iter
        E_old = E
        Ee = 0
        Fα = FermiMDzeros(Float64, (m,m))
        Fβ = FermiMDzeros(Float64, (m,m))
        H = T + V
        buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ints)
        F̃α = X*Fα*X
        F̃β = X*Fβ*X
        ϵα, C̃α = diagonalize(F̃α)
        ϵβ, C̃β = diagonalize(F̃β)
        Cα = X*C̃α
        Cβ = X*C̃β
        buildD!(Dα, Cα, Nα)
        buildD!(Dβ, Cβ, Nβ)
        for i in 1:m
            for j in 1:m
                Ee += 0.5 * (H[i,j]*(Dα[j,i]+Dβ[j,i]) + Fα[i,j]*Dα[j,i] + Fβ[i,j]*Dβ[j,i])
            end
        end
        E = Ee + molecule.Vnuc
        err = abs(E-E_old)
        println("E = $E,   ΔE = $err")
        ite += 1
        if abs(E - E_old) <= conv
            converged = true
            break
        end
    end
    ndocc = 4
    nsocc = 1
    nvir = 3
    ΔE = 0.00001
    Drms = 0.0002
    Orbitals = UHFOrbitals(molecule, ints.basis, ϵα, ϵβ, E, Cα, Cβ)
    return UHF(molecule, E, ndocc, nsocc, nvir, Orbitals, ΔE, Drms)
end

function buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ints)
    Fα .= H
    Fβ .= H
    ERI = ints["ERI"]
    calcJ!(Jα, Dα, ERI)
    calcJ!(Jβ, Dβ, ERI)
    calcK!(Kα, Dα, ERI)
    calcK!(Kβ, Dβ, ERI)
    Fα .+= Jα - Kα + Jβ
    Fβ .+= Jβ - Kβ + Jα
end

function calcJ!(J, D, ERI)
    @tensoropt J[i,j] = ERI[i,j,k,l] * D[l,k]
end

function calcK!(K, D, ERI)
    @tensoropt K[i,j] = ERI[i,l,k,j] * D[l,k]
end

function buildD!(D, C, N)
    Co = C[:,1:N]
    @tensoropt D[μ, ν] = Co[μ, i] * Co[ν, i]
end