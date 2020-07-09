module RHF
using Lints
using Fermi
using Fermi.Output
using LinearAlgebra

mutable struct RHFWavefunction{T} <: Fermi.AbstractHFWavefunction where T <: AbstractFloat
    refEnergy::T
    vnuc::T
    nocca::Int
    noccb::Int
    nvira::Int
    nvirb::Int
    basis::B where B <: Lints.BasisSet
    molecule::M where M <: Lints.Molecule
    Ca::Array{T,2} #AO->MO coefficients a
    Cb::Array{T,2} #AO->MO coefficients b
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    epsa::Array{T,1} #orbital eigenvalues a
    epsb::Array{T,1} #orbital eigenvalues b
    ERI::I where I <: Fermi.AbstractTensor#AO basis electron repulsion integrals
end

function RHFWavefunction(basis,molecule,nocca,noccb)
    RHFWavefunction(basis,molecule,nocca,noccb,
                          Fermi.ComputeEnvironment.interconnect,
                          Fermi.ComputeEnvironment.communicator,
                          Fermi.ComputeEnvironment.accelerator)
end

function RHFWavefunction(basis,molecule,nocca,noccb,
                               interconnect::Fermi.Environments.No_IC,
                               communicator::Fermi.Environments.NoCommunicator,
                               accelerator::Fermi.Environments.NoAccelerator)
    nprim = Lints.max_nprim(basis)
    l = Lints.max_l(basis)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,molecule)
    I_engines = []
    sz = Lints.getsize(S_engine,basis)
    for i in 1:Threads.nthreads()
        push!(I_engines,Lints.ERIEngine(nprim,l))
    end
    S = zeros(sz,sz)
    Ca = zeros(size(S))
    Cb = zeros(size(S))
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,basis)
    Lints.make_2D(T,T_engine,basis)
    Lints.make_2D(V,V_engine,basis)
    Lints.make_ERI(I,I_engines,basis)
    I = Fermi.MemTensor(I)
    ref = RHFWavefunction{Float64}(
                                   0.0,
                                   0.0,
                                   nocca,
                                   noccb,
                                   sz-nocca,
                                   sz-noccb,
                                   basis,
                                   molecule,
                                   Ca,
                                   Cb,
                                   S,
                                   T,
                                   V,
                                   zeros(Float64,sz),
                                   zeros(Float64,sz),
                                   I
                                  )
    RHFWavefunction(ref)

end

function RHFWavefunction(wfn::RHFWavefunction; doprint=false,maxit=50,Etol=1E-7,Dtol=1E-7)
    Fermi.HartreeFock.print_header()
    @output "    executing RHF\n"
    @output "    Forming initial Fock matrix ... "
    A = wfn.S^(-1/2)
    t = @elapsed begin
        Ft = transpose(A)*(wfn.T+wfn.V)*A
        e,Ct = eigen(Ft)
        C = A*Ct
        Co = C[:,1:wfn.nocca]
    end
    @output "done in {:>5.2f}s\n" t
    #G = 2*wfn.ERI.data - permutedims(wfn.ERI.data,(1,3,2,4))
    D = Fermi.contract(Co,Co,"um","vm")
    F = Fermi.contract(D,wfn.ERI,1.0,2.0,"rs","mnrs")
    Fermi.contract!(F,D,wfn.ERI,1.0,1.0,-1.0,"mn","rs","mrns")
    F += wfn.T
    F += wfn.V
    E = 0
    @output "\n"
    @output " Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "dE" "√|D|²" "t"
    @output repeat("~",80)*"\n"
    t = @elapsed for i in 1:maxit
        t_iter = @elapsed begin
            F .= 0
            F += wfn.T
            F += wfn.V
            Fermi.contract!(F,D,wfn.ERI,1.0,1.0,2.0,"mn","rs","mnrs")
            Fermi.contract!(F,D,wfn.ERI,1.0,1.0,-1.0,"mn","rs","mrns")
            Eelec = RHFEnergy(D,wfn.T+wfn.V,F)
            Enew = Eelec + wfn.vnuc
            Ft = transpose(A)*F*A
            Ft = Symmetric(Ft)
            e,Ct = eigen(Ft)#,sortby = x->-abs(x))
            C = A*Ct
            Co = C[:,1:wfn.nocca]
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
    wfn.epsa .= e
    wfn.epsb .= e
    wfn.Ca .= C 
    wfn.Cb .= C 
    wfn.refEnergy = E
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.17f}\n" E
    wfn
end
"""
    RHFEnergy
computes the RHF energy given a (not necessarily converged) core and
two electron integrals in the AO basis, as well as the desired density matrix.
## paramters
    h::Array{Float64,2} -> core hamiltonian array in AO basis
    g::Array{Float64,4} -> two electron integrals (eq 4b)
                        ->      g[μ,ν] = Σ(ρ,σ) ⟨μρ|νσ⟩ - ⟨μρ|σν⟩
    D::Array{Float64,2} -> density matrix (eq 2b) second C is actually Cconj
                        ->      D[μ,ν] = 2Σ(i;1:nocc) C[μ,i]*C[ν,i]
## outputs
    E::Float            -> SCF energy
"""
function RHFEnergy(D,H,F)
    sum(D .* (H .+ F))
end

function nmo(wfn::RHFWavefunction)
    return wfn.nocca + wfn.nvira
end
function Cao(wfn::RHFWavefunction)
    return wfn.Ca[:,1:wfn.nocca]
end
function Cav(wfn::RHFWavefunction)
    return wfn.Ca[:,wfn.nocca+1:wfn.nocca+wfn.nvira]
end
function Cbo(wfn::RHFWavefunction)
    return wfn.Cb[:,1:wfn.noccb]
end
function Cbv(wfn::RHFWavefunction)
    return wfn.Cb[:,wfn.noccb+1:wfn.noccb+wfn.nvirb]
end


end #module
