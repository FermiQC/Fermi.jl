function UHF(Alg::A) where A <: UHFAlgorithm
    ints = IntegralHelper{Float64}()
    UHF(ints, Alg)
end

function UHF(mol::Molecule, Alg::A) where A <: UHFAlgorithm
    UHF(IntegralHelper{Float64}(molecule=mol), Alg)
end

"""
    UHF(wfn::UHF, target_ints::IntegralHelper, Alg::A) where A <: UHFAlgorithm

Basis/geometry-projection guess (Werner 2004, same math as
`RHF(wfn::RHF, target_ints, Alg)` in RHFa.jl), applied separately to the alpha and
beta orbitals. Takes the target `IntegralHelper` explicitly so callers (e.g. the
geometry optimizer) can warm-start UHF SCF at a specific geometry without
mutating global options.
"""
function UHF(wfn::UHF, target_ints::IntegralHelper{Float64}, Alg::A) where A <: UHFAlgorithm

    Fermi.HartreeFock.uhf_header()
    output("Using {} wave function as initial guess", wfn.orbitals.basis)

    Sbb = target_ints["S"]
    Λ = Array(Sbb^(-1/2))

    bsA = GaussianBasis.BasisSet(wfn.orbitals.basis, wfn.molecule.atoms)
    bsB = GaussianBasis.BasisSet(target_ints.basis, target_ints.molecule.atoms)
    Sab = GaussianBasis.overlap(bsA, bsB)

    function project(Ca)
        T = transpose(Ca)*Sab*(Sbb^-1.0)*transpose(Sab)*Ca
        Cb = (Sbb^-1.0)*transpose(Sab)*Ca*T^(-1/2)
        return real.(Cb)
    end

    Cα = project(wfn.orbitals.Cα)
    Cβ = project(wfn.orbitals.Cβ)

    UHF(target_ints, Cα, Cβ, Λ, Alg)
end

function UHF(ints::IntegralHelper{Float64}, Alg::UHFa)
    Fermi.HartreeFock.uhf_header()
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
        Cα, Λ =  RHF_core_guess(ints)
        Cβ = deepcopy(Cα)
    elseif guess == "gwh"
        Cα, Λ = RHF_gwh_guess(ints)
        Cβ = deepcopy(Cα)
    end
    UHF(ints, Cα, Cβ, Λ, Alg)
end


function UHF(ints::IntegralHelper{Float64, <:AbstractERI, AtomicOrbitals}, Cα::AbstractMatrix, Cβ::AbstractMatrix, Λ::AbstractMatrix, Alg::UHFa)
    molecule = ints.molecule
    output(Fermi.string_repr(molecule))
    
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
    Vnuc = Molecules.nuclear_repulsion(molecule.atoms)
    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    ERI = ints["ERI"]
    m = size(S)[1]
    Dα = zeros(Float64, (m,m))
    Dβ = zeros(Float64, (m,m))
    Jα = zeros(Float64, (m,m))
    Jβ = zeros(Float64, (m,m))
    Kα = zeros(Float64, (m,m))
    Kβ = zeros(Float64, (m,m))
    ϵα = zeros(Float64, (m))
    ϵβ = zeros(Float64, (m))
    Fα = zeros(Float64, (m,m))
    Fβ = zeros(Float64, (m,m))
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


    build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints)
    output(" Guess Energy {:20.14f}", UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc))
 
    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[UHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    if do_diis
        DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end 
    t = @elapsed while ite <= maxit
        t_iter = @elapsed begin
            E_old = E
            
            # Transform Fock matrices to MO basis
            F̃α = Λ*Fα*Λ
            F̃β = Λ*Fβ*Λ
            
            # Solve for eigenvalues and eigenvectors
            ϵα, C̃α = LinearAlgebra.eigen(Symmetric(F̃α), sortby=x->x)
            ϵβ, C̃β = LinearAlgebra.eigen(Symmetric(F̃β), sortby=x->x)
            
            # Transform orbital coefficient matrices to AO basis
            Cα = Λ*C̃α
            Cβ = Λ*C̃β
            
            # Build density matrices
            buildD!(Dα, Cα, Nα)
            buildD!(Dβ, Cβ, Nβ)
            
            # Build Fock matrix
            build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints)
            
            # Calculate energy
            E = UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc)
            
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
    spinsquare=UHFSpin(Nα , Nβ, Cα, Cβ, S)

    output("    UHF done in {:>5.2f}s", t)
    output("    @Final UHF Energy     {:>20.12f} Eₕ", E)
    output("    <S^2> = {:> 15.10f}",spinsquare)
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
