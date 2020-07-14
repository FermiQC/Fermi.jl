function RHF()
    molecule = Molecule()
    RHF(molecule)
end

function RHF(Alg::ConventionalRHF)
    molecule = Molecule()
    aoint = ConventionalAOIntegrals(molecule)
    RHF(molecule, aoint, Alg)
end

function RHF(aoint::ConventionalAOIntegrals)
    molecule = Molecule()
    RHF(molecule, aoint)
end

function RHF(molecule::Molecule)
    aoint = ConventionalAOIntegrals(molecule)
    RHF(molecule, aoint)
end

function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals)
    if Fermi.CurrentOptions["scf_algorithm"] == "conventional"
        Alg = ConventionalRHF()
    else
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(Fermi.CurrentOptions["scf_algorithm"])"))
    end

    RHF(molecule, aoint, Alg)
end

function RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

    A = aoint.S^(-1/2)
    Ft = transpose(A)*(aoint.T+aoint.V)*A
    e,Ct = eigen(Ft)
    C = A*Ct

    RHF(molecule, aoint, C, Alg)
end

function RHF(wfn::RHF)

    if Fermi.CurrentOptions["scf_algorithm"] == "conventional"
        Alg = ConventionalRHF()
    else
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(Fermi.CurrentOptions["scf_algorithm"])"))
    end

    RHF(wfn, Alg)
end

function RHF(wfn::RHF, Alg::ConventionalRHF)

    aoint = ConventionalAOIntegrals(wfn.molecule)
    RHF(wfn, aoint, Alg)
end

function RHF(wfn::RHF, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

    RHF(wfn.molecule, aoint, wfn.C, Alg)
end

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
    @output " Nuclear Repulsion Energy:     {:5.10f}\n" molecule.Vnuc


    # Read in initial guess for orbitals
    x,y = size(Cguess)
    C = zeros(ndocc+nvir, ndocc+nvir)
    C[1:x, 1:y] .= Cguess

    # Form the orthogonalizer 
    A = aoint.S^(-1/2)

    # Form the density matrix
    Co = C[:, 1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)
    E = RHFEnergy(D, aoint.T+aoint.V, F) + molecule.Vnuc
    @output "\n"
    @output " Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t"
    @output repeat("~",80)*"\n"
    t = @elapsed for i in 1:maxit
        t_iter = @elapsed begin

            # Build the Fock Matrix
            build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)

            # Produce Ft
            Ft = transpose(A)*F*A

            # Get orbital energies and transformed coefficients
            e,Ct = eigen(Symmetric(Ft))

            # Reverse transformation to get MO coefficients
            C = A*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            Dnew = Fermi.contract(Co,Co,"um","vm")

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
        @output "    {:<3} {:>20.17f} {:>11.3e} {:>11.3e} {:>8.2f}\n" i E ΔE Drms t_iter

        if (ΔE < Etol) & (Drms < Dtol)
            build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)
            break
        end
    end
    @output repeat("~",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.10f}\n" E

    return RHF{Float64}(molecule, E, ndocc, nvir, C, F)
end

#function RHFEnergy(D::Array{Float64,2}, H::Array{Float64,2},F::Symmetric{Float64, Array{Float64,2}})
function RHFEnergy(D::Array{Float64,2}, H::Array{Float64,2},F::Array{Float64,2})
    return sum(D .* (H .+ F))
end

#function build_fock!(F::Symmetric{Float64, Array{Float64,2}}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Fermi.MemTensor)
function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Fermi.MemTensor)
    F .= H
    Fermi.contract!(F,D,ERI,1.0,1.0,2.0,"mn","rs","mnrs")
    Fermi.contract!(F,D,ERI,1.0,1.0,-1.0,"mn","rs","mrns")
end
