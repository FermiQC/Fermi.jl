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
        Cα, Λ =  UHF_core_guess(ints)
        Cβ = deepcopy(Cα)
    elseif guess == "gwh"
        Cα, Λ = UHF_gwh_guess(ints)
        Cβ = deepcopy(Cα)
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
    Dsα = deepcopy(Dα)
    Dsβ = deepcopy(Dβ)
    Fsα = deepcopy(Fα)
    Fsβ = deepcopy(Fβ)
    
    H = T + V
    Dα = buildD!(Dα, Cα, Nα)
    Dβ = buildD!(Dβ, Cβ, Nβ)
    Dα_old = deepcopy(Dα)
    Dβ_old = deepcopy(Dβ)
    Dsα = deepcopy(Dα)
    Dsβ = deepcopy(Dβ)


    buildfock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, H, Dα, Dβ, ERI)
    output(" Guess Energy {:20.14f}", UHFEnergy(H, Dα, Dβ, Fα, Fβ, molecule.Vnuc, m))
 
    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[UHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    if do_diis
        DMα = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        DMβ = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end 
    t = @elapsed while ite <= maxit
        t_iter = @elapsed begin
            E_old = E
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
            #Ee = 0
            #for i in 1:m
            #    for j in 1:m
            #        Ee += 0.5 * (H[i,j]*(Dα[j,i]+Dβ[j,i]) + Fα[i,j]*Dα[j,i] + Fβ[i,j]*Dβ[j,i])
            #    end
            #end
            #E = Ee + molecule.Vnuc
            E = UHFEnergy(H, Dα, Dβ, Fα, Fβ, molecule.Vnuc, m)
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

    nocc = Nα + Nβ  # TODO
    nvir = 2*m - nocc   #    output(repeat("-",80))

    output("    UHF done in {:>5.2f}s", t)
    output("    @Final UHF Energy     {:>20.12f} Eₕ", E)
    output("\n   • Orbitals Summary",)
    output("\n   ⬗ Alpha (α) orbitals")
    output("\n {:>10}   {:>15}   {:>10}", "Orbital", "Energy", "Occupancy")
    for i in eachindex(ϵα)
        output(" {:>10}   {:> 15.10f}   {:>6}", i, ϵα[i], (i ≤ Nα ? "↿" : ""))
    end
    output("\n   ⬗ Beta (β) orbitals")
    output("\n {:>10}   {:>15}   {:>10}", "Orbital", "Energy", "Occupancy")
    for i in eachindex(ϵβ)
        output(" {:>10}   {:> 15.10f}   {:>6}", i, ϵβ[i], (i ≤ Nβ ? "⇂" : ""))
    end
    output("")
    if converged
        output("   ✔  SCF Equations converged 😄")
    else
        output("❗ SCF Equations did not converge in {:>5} iterations ❗", maxit)
    end
    output(repeat("-",80))

    Orbitals = UHFOrbitals(molecule, ints.basis, ϵα, ϵβ, E, Cα, Cβ)
    return UHF(molecule, E, nocc, nvir, Orbitals, ΔE, Drms)
end
