function RMP2{T}(refWfn::Fermi.HartreeFock.RHF,alg::Fermi.MollerPlesset.Conventional) where T <: AbstractFloat
    print_header()
    dmp2 = 0.0
    nocc = refWfn.ndocc
    @output "*  executing restricted MP2\n"
    @output "   nocc: {:>3}\n" nocc
    @output "   nvir: {:>3}\n" refWfn.nvir
    rocc = 1:1:nocc
    rvir = nocc+1:1:nocc+refWfn.nvir
    epsa = refWfn.eps
    @output "*  performing AO->MO integral transformation ... "
    t = @elapsed moeri = refWfn.ints["OOVV"]
    @output "done in {:>5.2f}s\n" t
    @output "*  Computing MP2 energy ... "
    t = @elapsed for b in rvir
        for a in rvir
            for j in rocc
                for i in rocc
                    aa = a - nocc
                    bb = b - nocc
                    dmp2 +=
                        (
                            moeri[i, j, aa, bb] *
                            (2 * moeri[i, j, aa, bb] - moeri[i, j, bb, aa])
                        ) / (epsa[i] + epsa[j] - epsa[a] - epsa[b])
                end
            end
        end
    end
    @output "done in {:>5.2f}s\n" t
    @output "   @RMP2 {:>20.17f} Eh\n" dmp2
    RMP2(dmp2,Fermi.MemTensor{T}(zeros(0,0,0,0)))
end
