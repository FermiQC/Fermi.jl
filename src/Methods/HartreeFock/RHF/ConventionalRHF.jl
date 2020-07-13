function RHFWavefunction()
    molecule = Molecule()
    RHFWavefunction(molecule)
end

function RHFWavefunction(A::ConventionalRHF)
    molecule = Molecule()
    aoint = ConventionalAOIntegrals(molecule)
    RHFWavefunction(molecule, aoint, A)
end

function RHFWavefunction(aoint::ConventionalAOIntegrals)
    molecule = Molecule()
    RHFWavefunction(molecule, aoint)
end

function RHFWavefunction(molecule::Molecule)
    aoint = ConventionalAOIntegrals(molecule)
    RHFWavefunction(molecule, aoint)
end

function RHFWavefunction(molecule::Molecule, aoint::ConventionalAOIntegrals)
    if Fermi.CurrentOptions["scf_algorithm"] == "conventional"
        A = ConventionalRHF()
    else
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(Fermi.CurrentOptions["scf_algorithm"])"))
    end

    RHFWavefunction(molecule, aoint, A)
end

function RHFWavefunction(molecule::Molecule, aoint::ConventionalAOIntegrals, A::ConventionalRHF)
    Fermi.HartreeFock.print_header()

    maxit = Fermi.CurrentOptions["scf_max_iter"]
    Etol  = Fermi.CurrentOptions["e_conv"]
    Dtol  = Fermi.CurrentOptions["scf_max_rms"]

    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(Fermi.InvalidFermiOption("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end

    nvir = size(aoint.S)[1] - ndocc

    @output "    executing RHF\n"
    @output "    Forming initial Fock matrix ... "

    A = aoint.S^(-1/2)
    t = @elapsed begin
        Ft = transpose(A)*(aoint.T+aoint.V)*A
        e,Ct = eigen(Ft)
        C = A*Ct
        Co = C[:,1:ndocc]
    end
    @output "done in {:>5.2f}s\n" t
    #G = 2*wfn.ERI.data - permutedims(wfn.ERI.data,(1,3,2,4))
    D = Fermi.contract(Co,Co,"um","vm")
    F = Fermi.contract(D,aoint.ERI,1.0,2.0,"rs","mnrs")
    Fermi.contract!(F,D,aoint.ERI,1.0,1.0,-1.0,"mn","rs","mrns")
    F += aoint.T
    F += aoint.V
    E = 0
    @output "\n"
    @output " Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "dE" "√|D|²" "t"
    @output repeat("~",80)*"\n"
    t = @elapsed for i in 1:maxit
        t_iter = @elapsed begin
            F .= 0
            F += aoint.T
            F += aoint.V
            Fermi.contract!(F,D,aoint.ERI,1.0,1.0,2.0,"mn","rs","mnrs")
            Fermi.contract!(F,D,aoint.ERI,1.0,1.0,-1.0,"mn","rs","mrns")
            Eelec = RHFEnergy(D,aoint.T+aoint.V,F)
            Enew = Eelec + molecule.Vnuc
            Ft = transpose(A)*F*A
            Ft = Symmetric(Ft)
            e,Ct = eigen(Ft)#,sortby = x->-abs(x))
            C = A*Ct
            Co = C[:,1:ndocc]
            Dnew = Fermi.contract(Co,Co,"um","vm")
            dD = Dnew - D
            Drms = sqrt(sum(dD)^2)
            dE = Enew - E
            D = Dnew
            E = Enew
        end
        #if doprint println("@RHF $i $E $dE $Drms") end
        @output "    {:<3} {:>20.17f} {:>11.3e} {:>11.3e} {:>8.2f}\n" i E dE Drms t_iter
        if (dE < Etol) & (Drms < Dtol)
            break
        end
    end
    @output repeat("~",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.17f}\n" E

    return RHFWavefunction{Float64}(molecule, E, ndocc, nvir, C, e)
end

function RHFEnergy(D,H,F)
    sum(D .* (H .+ F))
end
