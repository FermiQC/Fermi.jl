function RMP2{T}(refwfn::RHF, ints::IntegralHelper{T}) where T <: AbstractFloat
    mp_header()

    # Check frozen core and inactive virtual
    nvir = refwfn.nvir
    ndocc = refwfn.ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(InvalidFermiOption("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(InvalidFermiOption("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end

    # Build range for MP2 amplitudes
    # Occupied range
    o = (1+core):ndocc
    # Virtual range
    v = (ndocc+1):(ndocc + nvir - inac)

    # Collect RHF information
    orbs = refwfn.orbitals
    C = T.(orbs.C)
    ϵo = T.(orbs.eps[o])
    ϵv = T.(orbs.eps[v])

    output("  Starting MP2 computation")
    output(" Number of frozen orbitals:             {:d}" , core)
    output(" Number of inactive orbitals:           {:d}" , inac)
    output(" Number of correlated electron pairs:   {:d}", length(o))
    output(" Number of correlated virtual orbitals: {:d}", length(v))
    output(" ⇒ Total number of MP2 amplitudes:      {:d}\n\n", length(o)^2*length(v)^2)

    output(repeat("-",80))

    if Options.get("mp2_type") == "conventional"
        output("• Performing AO->MO integral transformation...", ending="")
        t = @elapsed moeri = Fermi.Integrals.ao_to_mo_eri!(ints, C[:,o], C[:,v], C[:,o], C[:,v])
        output("   Done in {:>5.2f} s\n", t)

        output(" Computing MP2 Energy...", ending="")
        t = @elapsed begin
            Emp2 = zero(T)
            for b in 1:nvir
                for a in 1:nvir
                    for j in 1:ndocc
                        for i in 1:ndocc
                            Emp2 += moeri[i,a,j,b]*(2*moeri[i,a,j,b] - moeri[i,b,j,a]) /
                            (ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b])
                        end
                    end
                end
            end
        end
        output("   Done in {:>5.2f} s\n", t)
        output("   @Final RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
        output("   @Final RMP2 Total Energy       {:>20.12f} Eₕ", refwfn.energy+Emp2)
        output(repeat("-",80))

        RMP2(refwfn, Emp2, Emp2+refwfn.energy)

    elseif Options.get("mp2_type") == "df"
        
        output("• Performing AO->MO integral transformation...", ending="")
        t = @elapsed Bov = Fermi.Integrals.ao_to_mo_rieri!(ints, C[:,o], C[:,v])
        output("   Done in {:>5.2f} s\n", t)

        output(" Computing DF-MP2 Energy...", ending="")
        dfsz = size(Bov, 1)
        Bo_a = similar(Bov[:,:,1])
        Bo_b = similar(Bov[:,:,1])
        #t = @elapsed begin
        #    Emp2 = zero(T)
        #    for b in 1:nvir
        #        Bo_b .= Bov[:,:,b]
        #        for a in 1:nvir                        
        #            Bo_a .= Bov[:,:,a]
        #            @tensor moeri[i,j] := Bo_a[P,i]*Bo_b[P,j]
        #            for j in 1:ndocc
        #                for i in 1:ndocc
        #                    Emp2 += moeri[i,j]*(2*moeri[i,a,j,b] - moeri[i,b,j,a]) /
        #                    (ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b])
        #                end
        #            end
        #        end
        #    end
        #end
        output("   Done in {:>5.2f} s\n", t)
        output("   @Final DF-RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
        output("   @Final DF-RMP2 Total Energy       {:>20.12f} Eₕ", refwfn.energy+Emp2)
        output(repeat("-",80))

        RMP2(refwfn, Emp2, Emp2+refwfn.energy)
    end
end
