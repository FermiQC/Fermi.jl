"""
    do_rmp2

performs MP2 on a restricted HF reference.
## paramters
refWfn::Wfn -> wavefunction to which MP2 will be applied.

## outputs
dmp2::Float MP2 energy correction
"""
function do_rmp2(refWfn::Wfn;kwargs...)
    print_header()
    dmp2 = 0.0
    nocc = refWfn.nalpha
    @output "*  executing restricted MP2\n"
    @output "   nocc: {:>3}\n" refWfn.nalpha
    @output "   nvir: {:>3}\n" refWfn.nvira
    rocc = 1:1:refWfn.nalpha
    rvir = nocc+1:1:nocc+refWfn.nvira
    Cao = refWfn.Cao
    Cav = refWfn.Cav
    epsa = refWfn.epsa
    @output "*  performing AO->MO integral transformation ... "
    t = @elapsed moeri = get_eri(refWfn,"OOVV")#permutedims(tei_transform(refWfn.ao_eri, Cao, Cav, Cao, Cav, "oovv"), [1, 3, 2, 4])
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
    return dmp2
end
