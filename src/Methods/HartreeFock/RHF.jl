module RHF
using Fermi
using Fermi.Output
using LinearAlgebra

export do_RHF

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


function do_RHF(Fermi.ReferenceWavefunction; doprint=false,maxit=50,Etol=1E-7,Dtol=1E-7)
    Fermi.HartreeFock.print_header()
    @output "    executing RHF\n"
    #@output "    computing AO basis integrals ... \n"
    #@output "           using {:>2} engines for integral computation.\n" Threads.nthreads()
    #t = @elapsed begin 
    #    nprim = Lints.max_nprim(wfn.basis)
    #    l = Lints.max_l(wfn.basis)
    #    I_engines = []
    #    for i in 1:Threads.nthreads()
    #        push!(I_engines,Lints.ERIEngine(nprim,l))
    #    end
    #    Lints.make_ERI(wfn.ERI,I_engines,wfn.basis)
    #    G = Fermi.IntegralTransformation.antisymmetrize(wfn.ERI) #2*wfn.I - permutedims(wfn.I,[1,3,2,4])
    #end
    #@output "    done in {:>5.2f}s\n" t
    @output "    Forming initial Fock matrix ... "
    A = wfn.S^(-1/2)
    t = @elapsed begin
        Ft = transpose(A)*(wfn.T+wfn.V)*A
        e,Ct = eigen(Ft)
        C = wfn.A*Ct
        Co = C[:,1:wfn.ndocc]
    end
    @output "done in {:>5.2f}s\n" t
    #@tensor D[u,v] := Co[u,m]*Co[v,m]
    D = Fermi.contract(Co,Co,"um","vm")
    #@tensor F[m,n] := D[r,s]*G[m,n,r,s]
    F = Fermi.contract(D,G,"rs","mnrs")
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
            #@tensor F[m,n] := wfn.H[m,n] + D[r,s]*G[m,n,r,s]
            Fermi.contract!(F,D,G,"mn","rs","mnrs")
            Eelec = RHFEnergy(D,wfn.T+wfn.V,F)
            Enew = Eelec + wfn.vnuc
            Ft = transpose(A)*F*A
            Ft = Symmetric(Ft)
            e,Ct = eigen(Ft)#,sortby = x->-abs(x))
            C = A*Ct
            Co = C[:,1:wfn.nocca]
            #@tensor Dnew[u,v] := Co[u,m]*Co[v,m]
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

end #module
