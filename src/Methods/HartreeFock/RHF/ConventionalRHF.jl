using TensorOperations
"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function given Molecule, Integrals objects.
"""
function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

    @output "Using Core Guess\n"
    A = aoint.S^(-1/2)
    Ft = transpose(A)*(aoint.T+aoint.V)*A
    e,Ct = eigen(Ft)
    C = A*Ct

    RHF(molecule, aoint, C, Alg)
end

"""
    Fermi.HartreeFock.RHF(wfn::RHF, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function. Inital guess for orbitals is built from given RHF wfn.
"""
function RHF(wfn::RHF, Alg::ConventionalRHF)

    aoint = ConventionalAOIntegrals(wfn.molecule)
    RHF(wfn, aoint, Alg)
end

"""
    Fermi.HartreeFock.RHF(wfn::RHF, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function. Inital guess for orbitals is built from given RHF wfn. Integrals
are taken from the aoint input.
"""
function RHF(wfn::RHF, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

    # Projection of A→B done using equations described in Werner 2004 
    # https://doi.org/10.1080/0026897042000274801
    @output "Using {} wave function as initial guess\n" wfn.basis
    Ca = wfn.C
    Sbb = aoint.S
    Sab = Lints.projector(wfn.LintsBasis, aoint.LintsBasis)
    T = transpose(Ca)*Sab*(Sbb^-1)*transpose(Sab)*Ca
    Cb = (Sbb^-1)*transpose(Sab)*Ca*T^(-1/2)
    RHF(wfn.molecule, aoint, Cb, Alg)
end

"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Cguess::Array{Float64,2}, Alg::ConventionalRHF)

Basic function for conventional RHF using conventional integrals.
"""
function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, C::Array{Float64,2}, Alg::ConventionalRHF)
    # Print header
    do_diis = Fermi.CurrentOptions["diis"]
    Fermi.HartreeFock.print_header()

    do_diis ? DM = Fermi.DIIS.DIISManager{Float64}(size=8) : nothing 
    # Look up iteration options
    maxit = Fermi.CurrentOptions["scf_max_iter"]
    Etol  = 10.0^(-Fermi.CurrentOptions["e_conv"])
    Dtol  = Fermi.CurrentOptions["scf_max_rms"]

    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(Fermi.InvalidFermiOption("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end

    nvir = size(aoint.S)[1] - ndocc

    @output " Number of Doubly Occupied Orbitals:   {:5.0d}\n" ndocc
    @output " Number of Virtual Spatial Orbitals:   {:5.0d}\n" nvir

    # Form the orthogonalizer 
    A = aoint.S^(-1/2)

    # Form the density matrix
    Co = C[:, 1:ndocc]
    #D = Fermi.contract(Co,Co,"um","vm")
    @tensor D[u,v] := Co[u,m]*Co[v,m]
    
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)
    E = RHFEnergy(D, aoint.T+aoint.V, F) + molecule.Vnuc
    eps = Array{Float64, 1}(undef, ndocc+nvir)
    ite = 1
    converged = false

    display(aoint.S)
    @output "\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t" "DIIS"
    @output repeat("~",80)*"\n"

    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin

            # Build the Fock Matrix
            build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)

            # Produce Ft
            Ft = A*F*A
            #display(Ft)
            do_diis ? push!(DM, Ft) : nothing
            do_diis && ite > 3 ? Ft = Fermi.DIIS.extrapolate(DM) : nothing
            #display(Ft)

            # Get orbital energies and transformed coefficients
            eps,Ct = eigen((Ft))

            # Reverse transformation to get MO coefficients
            C = A*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            #Dnew = Fermi.contract(Co,Co,"um","vm")
            @tensor Dnew[u,v] := Co[u,m]*Co[v,m]

            # Compute Energy
            Eelec = RHFEnergy(Dnew, aoint.T+aoint.V, F)
            Enew = Eelec + molecule.Vnuc

            # Compute the Density RMS
            ΔD = Dnew - D
            Drms = sqrt(sum(ΔD.^2))

            # Compute Energy Change
            ΔE = Enew - E
            D .= Dnew
            E = Enew
        end
        @output "    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}\n" ite E ΔE Drms t_iter (do_diis && ite > 2)
        ite += 1

        if (abs(ΔE) < Etol) & (Drms < Dtol)
            converged = true
            break
        end
    end

    @output repeat("~",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.10f}\n" E

    @output "\n   • Orbitals Summary\n"
    @output "\n {:>10}   {:>15}   {:>10}\n" "Orbital" "Energy" "Occupancy"
    for i in eachindex(eps)
            @output " {:>10}   {:> 15.10f}   {:>6}\n" i eps[i] (i ≤ ndocc ? "↿⇂" : "")
    end
    
    if !converged
        @output "\n !! SCF Equations did not converge in {:>5} iterations !!\n" maxit
    end
    return RHF(aoint.basis, aoint.LintsBasis, molecule, E, ndocc, nvir, C, eps)
end

function RHFEnergy(D::Array{Float64,2}, H::Array{Float64,2},F::Array{Float64,2})
    return sum(D .* (H .+ F))
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Fermi.MemTensor)
    F .= H
    #Fermi.contract!(F,D,ERI,1.0,1.0,2.0,"mn","rs","mnrs")
    @tensoropt F[m,n] += 2.0*D[r,s]*ERI.data[m,n,r,s]
    #Fermi.contract!(F,D,ERI,1.0,1.0,-1.0,"mn","rs","mrns")
    @tensoropt F[m,n] -= D[r,s]*ERI.data[m,r,n,s]
end
