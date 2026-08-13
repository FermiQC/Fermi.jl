using GaussianBasis

function RHF(Alg::A) where A <: RHFAlgorithm
    warn_no_aoints()
    ints = IntegralHelper{Float64}()
    RHF(ints, Alg)
end

function RHF(mol::Molecule, Alg::A) where A <: RHFAlgorithm
    warn_no_aoints()
    RHF(IntegralHelper{Float64}(molecule=mol), Alg)
end

function RHF(ints::IntegralHelper{Float64}, Alg::A) where A <: RHFAlgorithm

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
        C, Λ = RHF_core_guess(ints)
    elseif guess == "gwh"
        C, Λ = RHF_gwh_guess(ints)
    end

    RHF(ints, C, Λ, Alg)
end

function RHF(wfn::RHF, Alg::A) where A <: RHFAlgorithm

    Fermi.HartreeFock.hf_header()

    # Projection of A→ B done using equations described in Werner 2004 
    # https://doi.org/10.1080/0026897042000274801

    output("Using {} wave function as initial guess", wfn.orbitals.basis)

    # B = target basis set
    warn_no_aoints()
    intsB = IntegralHelper{Float64}()

    # Assert both A and B have the same molecule.
    #if intsB.molecule != wfn.molecule
    #    output(" ! Input molecule does not match the molecule from the RHF wave function !")
    #end

    Sbb = intsB["S"]
    Λ = Array(Sbb^(-1/2))

    Ca = wfn.orbitals.C
    bsA = GaussianBasis.BasisSet(wfn.orbitals.basis, wfn.molecule.atoms)
    bsB = GaussianBasis.BasisSet(intsB.basis, intsB.molecule.atoms)
    Sab = GaussianBasis.overlap(bsA, bsB)

    T = transpose(Ca)*Sab*(Sbb^-1.0)*transpose(Sab)*Ca
    Cb = (Sbb^-1.0)*transpose(Sab)*Ca*T^(-1/2)
    Cb = real.(Cb)

    RHF(intsB, Cb, Λ, Alg)
end

"""
    RHF(wfn::RHF, target_ints::IntegralHelper, Alg::A) where A <: RHFAlgorithm

Same basis/geometry-projection guess as `RHF(wfn::RHF, Alg)` (Werner 2004), but
takes the target `IntegralHelper` explicitly instead of building one from the
current `molstring`/`basis` options -- useful for callers (e.g. the geometry
optimizer) that need to warm-start SCF at a specific, programmatically built
geometry without mutating global options every call.
"""
function RHF(wfn::RHF, target_ints::IntegralHelper{Float64}, Alg::A) where A <: RHFAlgorithm

    Fermi.HartreeFock.hf_header()
    output("Using {} wave function as initial guess", wfn.orbitals.basis)

    Sbb = target_ints["S"]
    Λ = Array(Sbb^(-1/2))

    Ca = wfn.orbitals.C
    bsA = GaussianBasis.BasisSet(wfn.orbitals.basis, wfn.molecule.atoms)
    bsB = GaussianBasis.BasisSet(target_ints.basis, target_ints.molecule.atoms)
    Sab = GaussianBasis.overlap(bsA, bsB)

    T = transpose(Ca)*Sab*(Sbb^-1.0)*transpose(Sab)*Ca
    Cb = (Sbb^-1.0)*transpose(Sab)*Ca*T^(-1/2)
    Cb = real.(Cb)

    RHF(target_ints, Cb, Λ, Alg)
end

function RHF(ints::IntegralHelper{Float64, <:AbstractERI, AtomicOrbitals}, C::AbstractMatrix, Λ::AbstractMatrix, Alg::RHFa)

    molecule = ints.molecule
    output(Fermi.string_repr(molecule))
    # Grab some options
    maxit = Options.get("scf_max_iter")
    Etol  = Options.get("scf_e_conv")
    Dtol  = Options.get("scf_max_rms")
    do_diis = Options.get("diis")
    oda = Options.get("oda")
    oda_cutoff = Options.get("oda_cutoff")
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

    # Grab ndocc,nvir,Vnuc
    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(FermiException("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end
    nvir = size(C,2) - ndocc
    nao = size(C,1)
    Vnuc = Molecules.nuclear_repulsion(molecule.atoms)

    output("Nuclear repulsion: {:15.10f}", Vnuc)
    output(" Number of AOs:                        {:5.0d}", nao)
    output(" Number of Doubly Occupied Orbitals:   {:5.0d}", ndocc)
    output(" Number of Virtual Spatial Orbitals:   {:5.0d}", nvir)

    
    S = ints["S"]
    T = ints["T"]
    V = ints["V"]
    ERI = ints["ERI"]

    # Form the density matrix from occupied subset of guess coeffs
    Co = C[:, 1:ndocc]
    @tensor D[u,v] := Co[u,m]*Co[v,m]
    D_old = deepcopy(D)
    eps = zeros(Float64,ndocc+nvir)

    # Build the inital Fock Matrix and diagonalize
    F = zeros(Float64,nao,nao)
    build_fock!(F, T + V, D, ints)
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
            eps, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)

            # Reverse transformation to get MO coefficients
            C = Λ*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            @tensor D[u,v] = Co[u,m]*Co[v,m]

            # Build the Fock Matrix
            build_fock!(F, T + V, D, ints)
            Eelec = RHFEnergy(D, T + V, F)

            # Compute Energy
            Enew = Eelec + Vnuc

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