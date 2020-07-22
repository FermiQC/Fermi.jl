using TensorOperations
"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function given Molecule, Integrals objects.
"""
function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

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

    RHF(wfn.molecule, aoint, wfn.C, Alg)
end

"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Cguess::Array{Float64,2}, Alg::ConventionalRHF)

Basic function for conventional RHF using conventional integrals.
"""
function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Cguess::Array{Float64,2}, Alg::ConventionalRHF)
    # Print header
    Fermi.HartreeFock.print_header()

    # Look up iteration options
    maxit = Fermi.CurrentOptions["scf_max_iter"]
    Etol  = Fermi.CurrentOptions["e_conv"]
    Dtol  = Fermi.CurrentOptions["scf_max_rms"]

    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(Fermi.InvalidFermiOption("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end

    nvir = size(aoint.S)[1] - ndocc

    @output " Number of Doubly Occupied Orbitals:   {:5.0d}\n" ndocc
    @output " Number of Virtual Spatial Orbitals:   {:5.0d}\n" nvir

    # Read in initial guess for orbitals
    x,y = size(Cguess)
    C = zeros(ndocc+nvir, ndocc+nvir)
    C[1:x, 1:y] .= Cguess

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

    @output "\n Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t"
    @output repeat("~",80)*"\n"

    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin

            # Build the Fock Matrix
            build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)

            # Produce Ft
            Ft = transpose(A)*F*A

            # Get orbital energies and transformed coefficients
            eps,Ct = eigen(Symmetric(Ft))

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
        @output "    {:<3} {:>20.17f} {:>11.3e} {:>11.3e} {:>8.2f}\n" ite E ΔE Drms t_iter
        ite += 1

        if (ΔE < Etol) & (Drms < Dtol)
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

    return RHF(molecule, E, ndocc, nvir, C, eps)
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
