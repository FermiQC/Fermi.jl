function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, Alg::A) where A <: RHFAlgorithm 
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
    nvir = size(aoint["S"])[1] - ndocc

    @output " Number of Doubly Occupied Orbitals:   {:5.0d}\n" ndocc
    @output " Number of Virtual Spatial Orbitals:   {:5.0d}\n" nvir
    
    # Form the orthogonalizer 
    sS = Hermitian(aoint["S"])
    Λ = sS^(-1/2)

    # Form the density matrix
    Co = C[:, 1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    F .= 0
    eps = Array{Float64, 1}(undef, ndocc+nvir)
    ite = 1
    converged = false
    T = aoint["T"]
    V = aoint["V"]
    Alg == ConventionalRHF() ? ERI = aoint["μ"] : nothing
    Alg == DFRHF()           ? ERI = aoint["B"] : nothing
    ΔE = 1.0
    Drms = 1.0
    oda_cutoff = 1E-1
    oda_shutoff = 20
    oda = Fermi.CurrentOptions["oda"]
    D_old = deepcopy(D)
    E = 0.0
    Drms = 1.0
    diis = false
    Co = C[:,1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")

    # Build the Fock Matrix
    build_fock!(F, T + V, D, ERI)
    Eelec = RHFEnergy(D, T + V, F)
    F̃ = deepcopy(F)
    D̃ = deepcopy(D)
    @output "\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t" "DIIS" "damp"
    @output repeat("~",80)*"\n"
    damp = 0.0 
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin



            #do_diis ? err = transpose(Λ)*(F*D*aoint["S"] - aoint["S"]*D*F)*Λ : nothing
            #do_diis ? push!(DM, F, err) : nothing
            #do_diis && ite > diis_start ? F = Fermi.DIIS.extrapolate(DM) : nothing


            # Produce Ft
            #if oda && Drms > oda_cutoff
            #    F = F̃
            #end
            if !oda || Drms < oda_cutoff
                Ft = Λ*F*transpose(Λ)
            else
                Ft = Λ*F̃*transpose(Λ)
            end

            # Get orbital energies and transformed coefficients

            eps,Ct = eigen(Hermitian(real.(Ft)))

            # Reverse transformation to get MO coefficients
            C = Λ*Ct
            Ct = real.(Ct)

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            D = Fermi.contract(Co,Co,"um","vm")

            # Build the Fock Matrix
            build_fock!(F, T + V, D, ERI)
            Eelec = RHFEnergy(D, T + V, F)
            #println(F ≈ F̃)
            #ite == 1 ? display(F) : nothing
            # Compute Energy
            Enew = Eelec + molecule.Vnuc


            if oda && Drms > oda_cutoff && ite < oda_shutoff
                diis = false
                dD = D̃ - D
                s = -tr(F̃ * dD)
                c = -tr((F - F̃) * (dD))
                if c <= -s/2
                    λ = 1.0
                else
                    λ = -s/(2*c)
                end
                F̃ .= (1-λ)*F̃ + λ*F
                D̃ .= (1-λ)*D̃ + λ*D
                damp = 1-λ
                do_diis ? err = transpose(Λ)*(F*D*aoint["S"] - aoint["S"]*D*F)*Λ : nothing
                do_diis ? push!(DM, F, err) : nothing
            end
            if (!oda || ite > oda_shutoff || Drms < oda_cutoff) && do_diis
                damp = 0.0
                diis = true
                D̃ = D
                do_diis ? err = transpose(Λ)*(F*D*aoint["S"] - aoint["S"]*D*F)*Λ : nothing
                do_diis ? push!(DM, F, err) : nothing
                do_diis && ite > diis_start ? F = Fermi.DIIS.extrapolate(DM) : nothing
            end
            #@tensor Dnew[u,v] := Co[u,m]*Co[v,m]

            # Compute the Density RMS
            ΔD = D - D_old
            Drms = sqrt(sum(ΔD.^2))

            # Compute Energy Change
            ΔE = Enew - E
            E = Enew
            D_old .= D
        end
        @output "    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.2f} {:>8}    {:5.4f}\n" ite E ΔE Drms t_iter diis damp
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

    aoint.C["O"] = C[:,1:ndocc]
    aoint.C["o"] = C[:,1:ndocc]
    aoint.C["V"] = C[:,ndocc+1:ndocc+nvir]
    aoint.C["v"] = C[:,ndocc+1:ndocc+nvir]
    aoint["F"]   = F

    return RHF(molecule, E, ndocc, nvir, eps, aoint)
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Array{Float64,4})
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
