function RHF(ints::IntegralHelper{Float64}, Alg::Direct)

    Fermi.HartreeFock.hf_header()

    output("Collecting One-electron integrals...")
    t = @elapsed begin
        ints["S"]
        ints["T"]
        ints["V"]
        if !haskey(ints.cache, "Jinv")
            ints["Jinv"] = inv(GaussianBasis.ERI_2e2c(ints.eri_type.basisset))
        end
    end
    output("Done in {:10.5f} s", t)

    guess = Options.get("scf_guess")
    if guess == "core"
        C, Λ = RHF_core_guess(ints)
    elseif guess == "gwh"
        C, Λ = RHF_gwh_guess(ints)
    end

    RHF(ints, C, Λ, Alg)
end

function RHF(ints::IntegralHelper{Float64, <:AbstractERI, AtomicOrbitals}, C::FermiMDArray{Float64,2}, Λ::FermiMDArray{Float64,2}, Alg::Direct)

    output("Using the DIRECT algorithm.")

    molecule = ints.molecule
    output(Fermi.string_repr(molecule))
    # Grab some options
    maxit       = Options.get("scf_max_iter")
    Etol        = Options.get("scf_e_conv")
    Dtol        = Options.get("scf_max_rms")
    do_diis     = Options.get("diis")
    oda         = Options.get("oda")
    oda_cutoff  = Options.get("oda_cutoff")
    oda_shutoff = Options.get("oda_shutoff")

    # Variables that will get updated iteration-to-iteration
    ite = 1
    E = 0.0
    ΔE = 1.0
    Drms = 1.0
    Drms = 1.0
    diis = false
    damp = 0.0 
    converged = false
    
    # Build a diis_manager, if needed
    if do_diis
        DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end

    # Grab ndocc,nvir
    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(FermiException("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end
    nvir = size(C,2) - ndocc
    nao = size(C,1)

    output(" Number of AOs:                        {:5.0d}", nao)
    output(" Number of Doubly Occupied Orbitals:   {:5.0d}", ndocc)
    output(" Number of Virtual Spatial Orbitals:   {:5.0d}", nvir)

    
    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    if !haskey(ints.cache, "Jinv")
        ints["Jinv"] = inv(GaussianBasis.ERI_2e2c(ints.eri_type.basisset))
    end
    Jinv = ints["Jinv"]

    # Form the density matrix from occupied subset of guess coeffs
    Co = C[:, 1:ndocc]
    D  = Co*Co'
    D_old = deepcopy(D)
    eps = nothing

    # Build the inital Fock Matrix
    F = T + V
    coulumb_to_fock!(F, D, Jinv, ints.orbitals.basisset, ints.eri_type.basisset)
    exchange_to_fock!(F, C, Jinv, ints.molecule.Nα, ints.orbitals.basisset, ints.eri_type.basisset)

    F̃ = deepcopy(F)
    D̃ = deepcopy(D)
    N = length(D) # Number of elements in D (For RMS computation)
    output(" Guess Energy {:20.14f}", RHFEnergy(D,T+V,F))

    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[RHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin
            # Produce Ft
            Ft = Λ'*F*Λ

            # Get orbital energies and transformed coefficients
            eps,Ct = diagonalize(Ft, hermitian=true)

            # Reverse transformation to get MO coefficients
            C = Λ*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            D = Co*Co'

            # Build the Fock Matrix
            F .= T + V
            coulumb_to_fock!(F, D, Jinv, ints.orbitals.basisset, ints.eri_type.basisset)
            exchange_to_fock!(F, C, Jinv, ints.molecule.Nα, ints.orbitals.basisset, ints.eri_type.basisset)
            Eelec = RHFEnergy(D, T + V, F)

            # Compute Energy
            Enew = Eelec + molecule.Vnuc

            # Store vectors for DIIS
            if do_diis
                err = transpose(Λ)*(F*D*S - S*D*F)*Λ
                push!(DM, F, err)
            end

            # Branch for ODA vs DIIS convergence aids
            diis = false
            damp = 0.0
            # Use ODA damping?
            if oda && Drms > oda_cutoff && ite < oda_shutoff
                diis = false
                dD = D - D̃
                s = tr(F̃ * dD)
                c = tr((F - F̃) * (dD))
                if c <= -s/(2*c)
                    λ = 1.0
                else
                    λ = -s/(2*c)
                end
                F̃ .= (1-λ)*F̃ + λ*F
                D̃ .= (1-λ)*D̃ + λ*D
                damp = 1-λ
                F .= F̃
            
            # Or Use DIIS?
            elseif do_diis && ite > diis_start
                diis = true
                F = Fermi.DIIS.extrapolate(DM)
            end

            # Compute the Density RMS
            ΔD = D - D_old
            Drms = √(sum(ΔD.^2) / N)

            # Compute Energy Change
            ΔE = Enew - E
            E = Enew
            D_old .= D
        end
        output("    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}    {:5.2f}", ite, E, ΔE, Drms, t_iter, diis, damp)
        ite += 1

        if (abs(ΔE) < Etol) & (Drms < Dtol) & (ite > 5)
            converged = true
            break
        end
    end

    output(repeat("-",80))
    output("    RHF done in {:>5.2f}s", t)
    output("    @Final RHF Energy     {:>20.12f} Eₕ", E)
    output("\n   • Orbitals Summary",)
    output("\n {:>10}   {:>15}   {:>10}", "Orbital", "Energy", "Occupancy")
    for i in eachindex(eps)
        output(" {:>10}   {:> 15.10f}   {:>6}", i, eps[i], (i ≤ ndocc ? "↿⇂" : ""))
    end
    output("")
    if converged
        output("   ✔  SCF Equations converged 😄")
    else
        output("❗ SCF Equations did not converge in {:>5} iterations ❗", maxit)
    end
    output(repeat("-",80))

    Orbitals = RHFOrbitals(molecule, ints.basis, eps, E, C)

    return RHF(molecule, E, ndocc, nvir, Orbitals, ΔE, Drms)
end