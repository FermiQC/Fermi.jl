function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, Alg::ConventionalRHF)
    RHF(molecule,aoint,C,aoint["μ"])
end
function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, Alg::DFRHF)
    RHF(molecule,aoint,C,aoint["B"])
end
"""
    RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64,N}) where N

The RHF kernel. Computes RHF on the given molecule with integral information defined in aoint. Starts from
the given C matrix. 
"""
function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64})
    Fermi.HartreeFock.print_header()

    #grab some options
    maxit = Fermi.CurrentOptions["scf_max_iter"]
    Etol  = 10.0^(-Fermi.CurrentOptions["e_conv"])
    Dtol  = Fermi.CurrentOptions["scf_max_rms"]
    do_diis = Fermi.CurrentOptions["diis"]
    oda = Fermi.CurrentOptions["oda"]
    oda_cutoff = Fermi.CurrentOptions["oda_cutoff"]
    oda_shutoff = Fermi.CurrentOptions["oda_shutoff"]

    #variables that will get updated iteration-to-iteration
    ite = 1
    E = 0.0
    ΔE = 1.0
    Drms = 1.0
    Drms = 1.0
    diis = false
    damp = 0.0 
    converged = false
    
    #build a diis_manager, if needed
    do_diis ? DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Fermi.CurrentOptions["ndiis"]) : nothing 
    do_diis ? diis_start = Fermi.CurrentOptions["diis_start"] : nothing

    #grab ndocc,nvir
    ndocc = try
        Int((molecule.Nα + molecule.Nβ)/2)
    catch InexactError
        throw(Fermi.InvalidFermiOption("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end
    nvir = size(aoint["S"])[1] - ndocc

    @output " Number of Doubly Occupied Orbitals:   {:5.0d}\n" ndocc
    @output " Number of Virtual Spatial Orbitals:   {:5.0d}\n" nvir
    
    # Form the orthogonalizer 
    S = Hermitian(aoint["S"])
    T = aoint["T"]
    V = aoint["V"]
    Λ = S^(-1/2)

    # Form the density matrix from occupied subset of guess coeffs
    Co = C[:, 1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    D_old = deepcopy(D)
    
    eps = zeros(Float64,ndocc+nvir)
    # Build the inital Fock Matrix and diagonalize
    F = zeros(Float64,ndocc+nvir,ndocc+nvir)
    build_fock!(F, T + V, D, ERI)
    F̃ = deepcopy(F)
    D̃ = deepcopy(D)

    @output "\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t" "DIIS" "damp"
    @output repeat("-",80)*"\n"
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin
            # Produce Ft
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

            # Compute Energy
            Enew = Eelec + molecule.Vnuc


            #branch for ODA vs DIIS convergence aids
            if oda && Drms > oda_cutoff && ite < oda_shutoff
                diis = false
                dD = D - D̃
                s = tr(F̃ * dD)
                c = tr((F - F̃) * (dD))
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
            elseif (!oda || ite > oda_shutoff || Drms < oda_cutoff) && do_diis
                damp = 0.0
                diis = true
                D̃ = D
                do_diis ? err = transpose(Λ)*(F*D*aoint["S"] - aoint["S"]*D*F)*Λ : nothing
                do_diis ? push!(DM, F, err) : nothing
                #do_diis && ite > diis_start ? F,_ = Fermi.DIIS.extrapolate(DM) : nothing
                if do_diis && ite > diis_start
                    F = Fermi.DIIS.extrapolate(DM)
                end
            end

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

    @output repeat("-",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.16f}\n" E
    @output "\n   • Orbitals Summary\n"
    @output "\n {:>10}   {:>15}   {:>10}\n" "Orbital" "Energy" "Occupancy"
    for i in eachindex(eps)
        @output " {:>10}   {:> 15.10f}   {:>6}\n" i eps[i] (i ≤ ndocc ? "↿⇂" : "")
    end
    @output "\n"
    if converged
        @output "   ✔  SCF Equations converged 😄\n" 
    else
        @output "❗ SCF Equations did not converge in {:>5} iterations ❗\n" maxit
    end
    @output repeat("-",80)*"\n"

    aoint.C["C"] = C
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
