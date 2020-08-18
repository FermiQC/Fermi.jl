
function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, Λ, Alg::ConventionalRHF)
    RHF(molecule,aoint,C,aoint["μ"],Λ)
end
function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, Λ, Alg::DFRHF)
    RHF(molecule,aoint,C,aoint["B"],Λ)
end
"""
    RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64,N}) where N

The RHF kernel. Computes RHF on the given molecule with integral information defined in aoint. Starts from
the given C matrix. 
"""
function RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64}, Λ::Array)
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
    
    S = Hermitian(aoint["S"])
    T = aoint["T"]
    V = aoint["V"]

    # Form the density matrix from occupied subset of guess coeffs
    Co = C[:, 1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    D_old = deepcopy(D)
    
    eps = zeros(Float64,ndocc+nvir)
    # Build the inital Fock Matrix and diagonalize
    F = zeros(Float64,ndocc+nvir,ndocc+nvir)
    build_fock!(F, T + V, D, ERI, Co)
    F̃ = deepcopy(F)
    D̃ = deepcopy(D)

    @output "\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}\n" "E[RHF]" "ΔE" "√|ΔD|²" "t" "DIIS" "damp"
    @output repeat("-",80)*"\n"
    t = @elapsed @fastmath while ite ≤ maxit
        t_iter = @elapsed begin
            # Produce Ft
            if !oda || Drms < oda_cutoff
                Ft = Λ'*F*Λ
            else
                Ft = Λ'*F̃*Λ
            end

            # Get orbital energies and transformed coefficients
            eps,Ct = eigen(Hermitian(real.(Ft)))

            # Reverse transformation to get MO coefficients
            C = Λ*Ct
            #Ct = real.(Ct)

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            #D = Fermi.contract(Co,Co,"um","vm")
            Fermi.contract!(D,Co,Co,0.0,1.0,1.0,"uv","um","vm")

            # Build the Fock Matrix
            build_fock!(F, T + V, D, ERI, Co)
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
            FPrms = sqrt(sum(err.^2))

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

    occ = CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in 1:ndocc])
    vir = CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in ndocc+1:ndocc+nvir])
    all = CanonicalOrbitals([CanonicalOrbital(Array{Float64,1}(C[:,i])) for i in 1:ndocc+nvir])
    aoint.orbs.ndocc = ndocc
    aoint.orbs.nvir = nvir
    aoint.orbs.frozencore = Fermi.CurrentOptions["drop_occ"]
    aoint.orbs.frozenvir = Fermi.CurrentOptions["drop_vir"]
    aoint.orbs["canonical"] = all
    aoint["F"]   = F

    return RHF(molecule, E, ndocc, nvir, eps, aoint)
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Array{Float64,4}, Co)
    F .= H
    Fermi.contract!(F,D,ERI,1.0,1.0,2.0,"mn","rs","mnrs")
    Fermi.contract!(F,D,ERI,1.0,1.0,-1.0,"mn","rs","mrns")
end

function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ERI::Array{Float64,3}, Co;build_K=true)
    F .= H
    sz = size(F,1)
    dfsz = size(ERI,1)
    no = size(Co,2)
    Fp = zeros(dfsz)
    Fermi.contract!(Fp,D,ERI,1.0,1.0,2.0,"Q","rs","Qrs")
    Fermi.contract!(F,Fp,ERI,1.0,1.0,1.0,"mn","Q","Qmn")
    ρ = zeros(dfsz,sz,no)
    Fermi.contract!(ρ,Co,ERI,"Qmc","rc","Qmr")
    Fermi.contract!(F,ρ,ρ,1.0,1.0,-1.0,"mn","Qmc","Qnc")
end

function build_J!(J::Array{Float64,2},D::Array{Float64,2},ERI::Array{Float64,3})
    sz = size(ERI,2)
    dfsz = size(ERI,1)
    Fp = zeros(dfsz)
    Fermi.contract!(Fp,D,ERI,1.0,1.0,2.0,"Q","rs","Qrs")
    Fermi.contract!(J,Fp,ERI,1.0,1.0,1.0,"mn","Q","Qmn")
end

function build_K!(K::Array{Float64,2},D::Array{Float64,2},ERI::Array{Float64,3})
    sz = size(ERI,2)
    dfsz = size(ERI,1)
    Fp = zeros(dfsz,sz,sz)
    Fermi.contract!(Fp,D,ERI,0.0,1.0,-1.0,"Qrn","rs","Qns")
    Fermi.contract!(K,ERI,Fp,1.0,1.0,1.0,"mn","Qmr","Qrn")
end
