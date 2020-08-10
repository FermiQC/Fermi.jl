function RMP2{T}(refwfn::Fermi.HartreeFock.RHF,alg::Fermi.MollerPlesset.DF) where T <: AbstractFloat
    Fermi.MollerPlesset.print_header()
    ints = refwfn.ints
    nocc    = refwfn.ndocc
    nvir    = refwfn.nvir
    eps = refwfn.eps
    ΔMP2 = 0.0
    @output "\tComputing MP2 with DF algorithm\n\n"
    @output "\tBasis: {}\n" ints.bname["primary"] 
    @output "\tDF basis: {}\n\n" ints.bname["aux"]
    ttotal = @elapsed begin
    @output "\tComputing and transforming integrals ... "
    t = @elapsed begin
        Fermi.Integrals.aux_ri!(ints)
        Bov = ints["BOV"]
    end
    @output " done in {:>5.2f} s\n" t
    Bi = zeros(T,size(Bov[:,1,:]))
    Bj = zeros(T,size(Bi))
    BAB = zeros(T,nvir,nvir)
    for i in 1:nocc
        for j in 1:nocc
            Bi[:,:] = Bov[:,i,:]
            Bj[:,:] = Bov[:,j,:]
            @tensoropt begin
                BAB[a,b] = Bi[Q,a]*Bj[Q,b]
            end
            BBA = transpose(BAB)
            for b in 1:nvir
                for a in 1:nvir
                    iajb = BAB[a,b]
                    ibja = BBA[a,b]
                    ΔMP2 += iajb*(2*iajb - ibja)/(eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc])
                end
            end
        end
    end
    #@tensor begin
    #    ΔMP2[] := Bov[Q,i,a]*Bov[Q,j,b]*(2*Bov[Q,i,a]*Bov[Q,j,b] - Bov[Q,i,b]*Bov[Q,j,a])/(eps[i] + eps[j] - eps[a] - eps[b])
    #end
    end
    @output "\n\tDF-MP2 energy is {:>14.10f}\n" ΔMP2
    @output "\tDF-MP2 done in {:>5.2f} s\n" ttotal
    @output repeat("-",80)*"\n"
    RMP2{T}(ΔMP2[],Fermi.MemTensor{T}(zeros(T,0,0,0,0)),Fermi.MemTensor{T}(zeros(T,0,0,0,0)))
end

