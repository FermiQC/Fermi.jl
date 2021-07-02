function UHF(Alg::UHFa)
    ints = IntegralHelper{Float64}(eri_type=Chonky())
    UHF(ints, Alg)
end

function UHF(ints::IntegralHelper{Float64}, Alg::UHFa)
    Fermi.HartreeFock.hf_header()
    output("Collecting necessary integrals...")
    t = @elapsed begin
        ints["S"]
        ints["T"]
        ints["V"]
        ints["ERI"]
    end
    output("Done in {:10.5f} s", t)

    guess = Options.get("scf_guess")
    if guess == "core"
        throw(FermiException("Not for UHF bub"))
    elseif guess == "gwh"
        throw(FermiException("Not for UHF bub"))
    elseif guess == "no"
        S = ints["S"]
        m = size(S)[1]
        Cα = FermiMDzeros(Float64, (m,m))
        Cβ = FermiMDzeros(Float64, (m,m))
        Λ = S^(-1/2)
    end
    UHF(ints, Cα, Cβ, Λ, Alg)
end


function UHF(ints::IntegralHelper{Float64, <:AbstractERI, AtomicOrbitals}, Cα::FermiMDArray{Float64, 2}, Cβ::FermiMDArray{Float64, 2}, Λ::FermiMDArray{Float64, 2}, Alg::UHFa)
    molecule = ints.molecule
    output(Fermi.Geometry.string_repr(molecule))
    
    # Grab options
    maxit = Options.get("scf_max_iter")
    Etol = Options.get("scf_e_conv")
    Dtol = Options.get("scf_max_rms")
    do_diis = Options.get("diis")
    oda = Options.get("oda")
    oda_cutoff = Options.get("oda_cutoff")
    oda_shutoff = Options.get("oda_shutoff")

    # Loop variables
    ite = 1
    E = 0.0
    ΔE = 1.0
    Drms = 1.0
    diis = false
    damp = 0.0
    converged = false
    
    Nα = molecule.Nα
    Nβ = molecule.Nβ
    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    ERI = ints["ERI"]
    m = size(S)[1]
    Dα = FermiMDzeros(Float64, (m,m))
    Dβ = FermiMDzeros(Float64, (m,m))
    Jα = FermiMDzeros(Float64, (m,m))
    Jβ = FermiMDzeros(Float64, (m,m))
    Kα = FermiMDzeros(Float64, (m,m))
    Kβ = FermiMDzeros(Float64, (m,m))
    ϵα = FermiMDzeros(Float64, (m))
    ϵβ = FermiMDzeros(Float64, (m))
    Fα = FermiMDzeros(Float64, (m,m))
    Fβ = FermiMDzeros(Float64, (m,m))
    Dα_old = deepcopy(Dα)
    Dβ_old = deepcopy(Dβ)
    Dsα = deepcopy(Dα)
    Dsβ = deepcopy(Dβ)
    Fsα = deepcopy(Fα)
    Fsβ = deepcopy(Fβ)
    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[RHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    if do_diis
        DMα = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        DMβ = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end 
    t = @elapsed while ite <= maxit
        t_iter = @elapsed begin
            E_old = E
            H = T + V
            buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ERI)
            # Transform Fock matrices to MO basis
            F̃α = Λ*Fα*Λ
            F̃β = Λ*Fβ*Λ
            # Solve for eigenvalues and eigenvectors
            ϵα, C̃α = diagonalize(F̃α)
            ϵβ, C̃β = diagonalize(F̃β)
            # Transform orbital coefficient matrices to AO basis
            Cα = Λ*C̃α
            Cβ = Λ*C̃β
            # Build density matrices
            buildD!(Dα, Cα, Nα)
            buildD!(Dβ, Cβ, Nβ)
            # Calculate energy
            Ee = 0
            for i in 1:m
                for j in 1:m
                    Ee += 0.5 * (H[i,j]*(Dα[j,i]+Dβ[j,i]) + Fα[i,j]*Dα[j,i] + Fβ[i,j]*Dβ[j,i])
                end
            end
            E = Ee + molecule.Vnuc
            # Store vectors for DIIS
            if do_diis
                err_α = transpose(Λ)*(Fα*Dα*S - S*Dα*Fα)*Λ
                err_β = transpose(Λ)*(Fβ*Dβ*S - S*Dβ*Fβ)*Λ
                push!(DMα, Fα, err_α)
                push!(DMβ, Fβ, err_β)
            end

            # Branch for ODA vs DIIS convergence aids
            diis = false
            damp = 0.0
            # Use ODA damping?
            if oda && Drms > oda_cutoff && ite < oda_shutoff
                diis = false
                odadamping!(diis, damp, Dα, Dsα, Fα, Fsα)
                odadamping!(diis, damp, Dβ, Dsβ, Fβ, Fsβ)
            # Or Use DIIS?
            elseif do_diis && ite > diis_start
                diis = true
                Fα = Fermi.DIIS.extrapolate(DMα)
                Fβ = Fermi.DIIS.extrapolate(DMβ)
            end

            # Calculate energy difference, Drms, and check for convergence
            ΔE = E-E_old
            ΔDα = Dα - Dα_old 
            ΔDβ = Dβ - Dβ_old
            Drms = (sum(ΔDα.^2)/m^2)^(1/2) + (sum(ΔDβ.^2)/m^2)^(1/2)
            Dα_old .= Dα
            Dβ_old .= Dβ
        end
        output("    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}    {:5.2f}", ite, E, ΔE, Drms, t_iter, diis, damp)
        ite += 1
        if (abs(ΔE) <= Etol) & (Drms <= Dtol) & (ite > 5)
            converged = true
            break
        end
    end
    ndocc = 4  # TODO
    nsocc = 1  #
    nvir = 3   #
    Orbitals = UHFOrbitals(molecule, ints.basis, ϵα, ϵβ, E, Cα, Cβ)
    return UHF(molecule, E, ndocc, nsocc, nvir, Orbitals, ΔE, Drms)
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