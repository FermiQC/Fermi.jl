function RHF(molecule::Molecule, aoint::I, C::Array{Float64,2}, Alg::A) where { I <: AbstractAOIntegrals,
                                                                                A <: RHFAlgorithm }
    Fermi.HartreeFock.print_header()
    do_diis = Fermi.CurrentOptions["diis"]
    do_diis ? DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Fermi.CurrentOptions["ndiis"]) : nothing 
    do_diis ? diis_start = Fermi.CurrentOptions["diis_start"] : nothing
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
    sS = Hermitian(aoint.S)
    Λ = sS^(-1/2)

    # Form the density matrix
    Co = C[:, 1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    eps = Array{Float64, 1}(undef, ndocc+nvir)
    ite = 1
    converged = false

    E = 0.0
    Drms = 1.0
    @output "\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t" "DIIS"
    @output repeat("~",80)*"\n"
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin

            # Build the Fock Matrix
            F_old = deepcopy(F)
            build_fock!(F, aoint.T+aoint.V, D, aoint.ERI)
            Eelec = RHFEnergy(D, aoint.T+aoint.V, F)
            #ite == 1 ? display(F) : nothing


            do_diis ? err = transpose(Λ)*(F*D*aoint.S - aoint.S*D*F)*Λ : nothing
            do_diis ? push!(DM, F, err) : nothing
            do_diis && ite > diis_start ? F = Fermi.DIIS.extrapolate(DM) : nothing

            # Produce Ft
            Ft = Λ*F*transpose(Λ)

            # Get orbital energies and transformed coefficients
            eps,Ct = eigen(Hermitian(Ft))

            # Reverse transformation to get MO coefficients
            C = Λ*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            Dnew = Fermi.contract(Co,Co,"um","vm")
            #@tensor Dnew[u,v] := Co[u,m]*Co[v,m]

            # Compute Energy
            Enew = Eelec + molecule.Vnuc

            # Compute the Density RMS
            ΔD = Dnew - D
            Drms = sqrt(sum(ΔD.^2))

            # Compute Energy Change
            ΔE = Enew - E
            E = Enew
            D .= Dnew 
        end
        @output "    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}\n" ite E ΔE Drms t_iter (do_diis && ite > diis_start)
        ite += 1

        if (abs(ΔE) < Etol) & (Drms < Dtol) & (ite > 5)
            converged = true
            break
        end
    end

    @output repeat("~",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.16f}\n" E

    @output "\n   • Orbitals Summary\n"
    @output "\n {:>10}   {:>15}   {:>10}\n" "Orbital" "Energy" "Occupancy"
    for i in eachindex(eps)
            @output " {:>10}   {:> 15.10f}   {:>6}\n" i eps[i] (i ≤ ndocc ? "↿⇂" : "")
    end
    
    if !converged
        @output "\n !! SCF Equations did not converge in {:>5} iterations !!\n" maxit
    end
    return RHF(aoint.bname, molecule, E, ndocc, nvir, C, eps, aoint)
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Fermi.MemTensor)
    F .= H
    Fermi.contract!(F,D,ERI,1.0,1.0,2.0,"mn","rs","mnrs")
    Fermi.contract!(F,D,ERI,1.0,1.0,-1.0,"mn","rs","mrns")
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Array{Float64,3})
    F .= H
    sz = size(F,1)
    dfsz = size(ERI,1)
    Fp = zeros(dfsz)
    Fermi.contract!(Fp,D,ERI,1.0,1.0,2.0,"Q","rs","Qrs")
    Fermi.contract!(F,Fp,ERI,1.0,1.0,1.0,"mn","Q","Qmn")
    Fp = zeros(dfsz,sz,sz)
    Fermi.contract!(Fp,D,ERI,0.0,1.0,-1.0,"Qrn","rs","Qns")
    Fermi.contract!(F,ERI,Fp,1.0,1.0,1.0,"mn","Qmr","Qrn")
end
