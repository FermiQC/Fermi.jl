function RMP2{T}(refwfn::RHF, ints::IntegralHelper{T}) where T <: AbstractFloat
    mp_header()

    orbs = refwfn.orbitals
    C = T.(orbs.C)
    # Occupied range
    o = 1:refwfn.ndocc
    # Virtual range
    v = (refwfn.ndocc+1):(refwfn.ndocc+refwfn.nvir)

    ϵo = T.(orbs.eps[o])
    ϵv = T.(orbs.eps[v])

    output("*  executing restricted MP2")
    output("   nocc: {:>3}\n" ,refwfn.ndocc)
    output("   nvir: {:>3}\n" ,refwfn.nvir)

    output("*  performing AO->MO integral transformation... ")
    t = @elapsed moeri = Fermi.Integrals.ao_to_mo_eri!(ints, C[:,o], C[:,v], C[:,o], C[:,v])
    output("done in {:>5.2f}s\n", t)

    mp2type = "greedy"

    if mp2type == "greedy"
        D = [1/(i+j-a-b) for i = ϵo, a = ϵv, j = ϵo, b = ϵv]

        V = 2*moeri - permutedims(moeri, (1,4,3,2))

        #@tensoropt Emp2 = moeri[i,a,j,b]*V[i,a,j,b]*D[i,a,j,b]
        Emp2 = sum(moeri.*V.*D)
        output("   @RMP2 {:>20.17f} Eh\n", Emp2)
        output("   @RMP2 {:>20.17f} Eh\n", refwfn.energy+Emp2)
        println(typeof(Emp2))
    end

    #@output "*  Computing MP2 energy ... "
    #t = @elapsed for b in rvir
    #    for a in rvir
    #        for j in rocc
    #            for i in rocc
    #                aa = a - nocc
    #                bb = b - nocc
    #                dmp2 +=
    #                    (
    #                        moeri[i, j, aa, bb] *
    #                        (2 * moeri[i, j, aa, bb] - moeri[i, j, bb, aa])
    #                    ) / (epsa[i] + epsa[j] - epsa[a] - epsa[b])
    #            end
    #        end
    #    end
    #end
    #@output "done in {:>5.2f}s\n" t
    #RMP2(dmp2,Fermi.MemTensor{T}(zeros(0,0,0,0)))
end
