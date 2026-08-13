function UHF(ints::IntegralHelper{Float64}, Alg::UHFDirect)

    Fermi.HartreeFock.uhf_header()

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
        Cα, Λ = RHF_core_guess(ints)
    elseif guess == "gwh"
        Cα, Λ = RHF_gwh_guess(ints)
    end
    Cβ = deepcopy(Cα)

    UHF(ints, Cα, Cβ, Λ, Alg)
end

function UHF(ints::IntegralHelper{Float64, <:AbstractERI, AtomicOrbitals}, Cα::AbstractMatrix, Cβ::AbstractMatrix, Λ::AbstractMatrix, Alg::UHFDirect)

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
    diis = false
    damp = 0.0
    converged = false

    # Build a diis_manager, if needed
    if do_diis
        DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end

    Nα = molecule.Nα
    Nβ = molecule.Nβ
    Vnuc = Molecules.nuclear_repulsion(molecule.atoms)

    output("Nuclear repulsion: {:15.10f}", Vnuc)

    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    if !haskey(ints.cache, "Jinv")
        ints["Jinv"] = inv(GaussianBasis.ERI_2e2c(ints.eri_type.basisset))
    end
    Jinv = ints["Jinv"]
    H = T + V
    m = size(S,1)

    # Build the density matrices from the occupied subset of guess coefficients
    Dα = zeros(Float64, m, m)
    Dβ = zeros(Float64, m, m)
    buildD!(Dα, Cα, Nα)
    buildD!(Dβ, Cβ, Nβ)
    Dα_old = deepcopy(Dα)
    Dβ_old = deepcopy(Dβ)
    Dsα = deepcopy(Dα)
    Dsβ = deepcopy(Dβ)

    Fα = zeros(Float64, m, m)
    Fβ = zeros(Float64, m, m)
    Fsα = deepcopy(Fα)
    Fsβ = deepcopy(Fβ)

    b1 = ints.orbitals.basisset
    b2 = ints.eri_type.basisset

    # Coulomb is linear in the density, so it's built once from the total density
    # (Dα+Dβ) and shared by both spins, instead of computing J[Dα] and J[Dβ]
    # separately as the in-core/DF build_fock! does. Exchange stays spin-separated
    # since it depends on the same-spin occupied orbitals only.
    function build_fock_direct!(Fα, Fβ, Dα, Dβ, Cα, Cβ)
        Fα .= H
        # coulumb_to_fock! expects a "half density" (the RHF convention D = Co*Co',
        # doubled internally by its own 2.0 factor to reach the true physical
        # density). Dα+Dβ is already the true total physical density for UHF, so
        # it must be halved here for that internal doubling to land correctly.
        coulumb_to_fock!(Fα, (Dα .+ Dβ) ./ 2, Jinv, b1, b2)
        Fβ .= Fα
        exchange_to_fock!(Fα, Cα, Jinv, Nα, b1, b2)
        exchange_to_fock!(Fβ, Cβ, Jinv, Nβ, b1, b2)
    end

    build_fock_direct!(Fα, Fβ, Dα, Dβ, Cα, Cβ)
    output(" Guess Energy {:20.14f}", UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc))

    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[UHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))

    ϵα = nothing
    ϵβ = nothing
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin
            E_old = E

            # Transform Fock matrices to MO basis
            F̃α = Λ'*Fα*Λ
            F̃β = Λ'*Fβ*Λ

            # Get orbital energies and transformed coefficients
            ϵα, C̃α = LinearAlgebra.eigen(Symmetric(F̃α), sortby=x->x)
            ϵβ, C̃β = LinearAlgebra.eigen(Symmetric(F̃β), sortby=x->x)

            # Reverse transformation to get MO coefficients
            Cα = Λ*C̃α
            Cβ = Λ*C̃β

            # Produce new density matrices
            buildD!(Dα, Cα, Nα)
            buildD!(Dβ, Cβ, Nβ)

            # Build the Fock matrices
            build_fock_direct!(Fα, Fβ, Dα, Dβ, Cα, Cβ)
            Enew = UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc)

            # Store vectors for DIIS
            if do_diis
                err_α = transpose(Λ)*(Fα*Dα*S - S*Dα*Fα)*Λ
                err_β = transpose(Λ)*(Fβ*Dβ*S - S*Dβ*Fβ)*Λ
                err_v = vcat(err_α, err_β)
                F_v = vcat(Fα, Fβ)
                push!(DM, F_v, err_v)
            end

            # Branch for ODA vs DIIS convergence aids
            diis = false
            damp = 0.0
            # Use ODA damping?
            if oda && Drms > oda_cutoff && ite < oda_shutoff
                damp = odadamping(Dα, Dsα, Fα, Fsα)
                damp = odadamping(Dβ, Dsβ, Fβ, Fsβ)
            # Or Use DIIS?
            elseif do_diis && ite > diis_start
                diis = true
                F_v = Fermi.DIIS.extrapolate(DM)
                Fα .= F_v[1:m, :]
                Fβ .= F_v[m+1:2m, :]
            end

            # Compute the Density RMS
            ΔDα = Dα - Dα_old
            ΔDβ = Dβ - Dβ_old
            Drms = (sum(ΔDα.^2)/m^2)^(1/2) + (sum(ΔDβ.^2)/m^2)^(1/2)

            # Compute Energy Change
            ΔE = Enew - E
            E = Enew
            Dα_old .= Dα
            Dβ_old .= Dβ
        end
        output("    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}    {:5.2f}", ite, E, ΔE, Drms, t_iter, diis, damp)
        ite += 1

        if (abs(ΔE) < Etol) & (Drms < Dtol) & (ite > 5)
            converged = true
            break
        end
    end

    nocc = Nα + Nβ
    nvir = 2*m - nocc

    output(repeat("-",80))
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
