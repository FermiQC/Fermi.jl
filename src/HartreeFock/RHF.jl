module RHF
using Fermi
using TensorOperations
using LinearAlgebra
using Lints
using Fermi.Output

export RHF
export RHFWfn
export RHFCompute

struct RHFWfn
    energy::Array{Float64,1}
    molecule::Lints.MoleculeAllocated
    basis::Lints.BasisSetAllocated
    nelec::Int
    C::Array{Float64,2}
    eps::Array{Float64,1}
    H::Array{Float64,2}
    S::Array{Float64,2}
    A::Array{Float64,2}
    D::Array{Float64,2}
    I::Array{Float64,4}
    vnuc::Float64
    ndocc::Int
    grad::Bool
    hess::Bool
    GradN::Array{Float64,2}
    GradS::Array{Float64,2}
    GradSp::Array{Float64,2}
    GradV::Array{Float64,2}
    GradT::Array{Float64,2}
    GradJ::Array{Float64,2}
    GradK::Array{Float64,2}
    Grad::Array{Float64,2}
end

function RHFWfn(basis,molecule,nelec;debug=false,grad=false,hess=false)
    dummy2 = Array{Float64}(undef,0,0)
    dummy4 = Array{Float64}(undef,0,0,0,0)
    vnuc = 0
    get_nuclear_repulsion()
    nprim = Lints.max_nprim(basis)
    l = Lints.max_l(basis)
    E = zeros(Float64,1)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,molecule)
    I_engine = Lints.ERIEngine(nprim,l)
    sz = Lints.getsize(S_engine,basis)
    S = zeros(sz,sz)
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,basis)
    Lints.make_2D(T,T_engine,basis)
    Lints.make_2D(V,V_engine,basis)
    A = S^(-1/2)
    H = T+V
    C = zeros(sz,sz)
    D = zeros(sz,sz)
    if ! grad
        GradN = dummy2
        GradS = dummy2
        GradSp = dummy2
        GradV = dummy2
        GradT = dummy2
        GradJ = dummy2
        GradK = dummy2
        Grad = dummy2
    end
    RHFWfn(E,molecule,basis,nelec,C,zeros(sz),H,S,A,D,I,vnuc,nelec/2,grad,hess,GradN,GradS,
          GradSp,GradV,GradT,GradJ,GradK,Grad)
end


function RHFCompute(wfn::RHFWfn;doprint=false,maxit=50,Etol=1E-7,Dtol=1E-7)
    Fermi.HartreeFock.print_header()
    @output "    executing RHF\n"
    @output "    computing AO basis integrals ... \n"
    @output "           using {:>2} engines for integral computation.\n" Threads.nthreads()
    nprim = Lints.max_nprim(wfn.basis)
    l = Lints.max_l(wfn.basis)
    engines = []
    for i in 1:Threads.nthreads()
        push!(engines,Lints.ERIEngine(nprim,l))
    end
    Lints.make_ERI(wfn.I,engines,wfn.basis)
    G = 2*wfn.I - permutedims(wfn.I,[1,3,2,4])
    t = 0.0
    @output "    done in {:>5.2f}s\n" t
    @output "    Forming initial Fock matrix ... "
    t = @elapsed begin
        Ft = transpose(wfn.A)*wfn.H*wfn.A
        e,Ct = eigen(Ft)
        C = wfn.A*Ct
        Co = C[:,1:wfn.ndocc]
    end
    @output "done in {:>5.2f}s\n" t
    @tensor D[u,v] := Co[u,m]*Co[v,m]
    @tensor F[m,n] := D[r,s]*G[m,n,r,s]
    F += wfn.H
    E = 0#RHFEnergy(D,wfn.H,F) + wfn.vnuc
    #if doprint println("@RHF 0 $E") end

    @output "\n"
    @output " Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "dE" "√|D|²" "t"
    @output repeat("~",80)*"\n"
    t = @elapsed for i in 1:maxit
        t_iter = @elapsed begin
            @tensor F[m,n] := wfn.H[m,n] + D[r,s]*G[m,n,r,s]
            Eelec = RHFEnergy(D,wfn.H,F)
            Enew = Eelec + wfn.vnuc
            Ft = transpose(wfn.A)*F*wfn.A
            Ft = Symmetric(Ft)
            e,Ct = eigen(Ft)#,sortby = x->-abs(x))
            C = wfn.A*Ct
            Co = C[:,1:wfn.ndocc]
            @tensor Dnew[u,v] := Co[u,m]*Co[v,m]
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
    wfn.eps .= e
    wfn.C .= C 
    wfn.energy[1] = E
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.17f}\n" E
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
    temp = H + F
    @tensor begin
        E[] := D[m,n]*temp[m,n]
    end
    return E[]
end

function get_nuclear_repulsion()
    open("/tmp/molfile.xyz","w") do molfile
        for l in eachline(molfile)
            println(l)
        end
    end
end

end #module
