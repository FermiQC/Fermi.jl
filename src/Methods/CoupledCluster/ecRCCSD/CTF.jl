"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function ecRCCSD{T}(Alg::CTF) where T <: AbstractFloat
    @output "Calling CASCI module...\n"
    # Call CASCI
    cas = Fermi.ConfigurationInteraction.CASCI{T}()

    # Print intro
    Fermi.CoupledCluster.print_header()
    @output "\n    • Computing Externally Corrected CCSD with the ecRCCSD module.\n\n"

    # Save reference wavefunction and process CAS data. Modify Ref (if not HF)
    @output "Processing CAS data...\n"
    Casdata = get_cas_data(cas)
    refwfn = cas.ref
    ints = refwfn.ints
    # Delete data that will not be used
    cas = nothing

    # Get MO Integrals
    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    # Process CAS data to get ecT1 and ecT2 (Cluster Decomposition step)
    frozen = Fermi.CurrentOptions["cas_frozen"]
    active = Fermi.CurrentOptions["cas_active"] ≢ -1 ? Fermi.CurrentOptions["cas_active"] : refwfn.nvir+refwfn.ndocc-frozen

    if drop_occ > frozen
        error("\nFrozen orbitals in the CC step ($drop_occ) cannot be greater than the number of frozen orbitals in the CAS ($frozen).")
    end

    if drop_vir > (refwfn.ndocc+refwfn.nvir) - active - frozen
        error("\nToo many virtual orbitals dropped ($drop_vir) for the active space.")
    end

    actocc = collect((1+frozen):refwfn.ndocc)
    actvir = collect((1+refwfn.ndocc):(active+frozen))

    @output "\n   • CAS Decomposition started:\n"
    @output "Active Occupied Orbitals: {}\n" actocc
    @output "Active Virtual Orbitals:  {}\n" actvir
    t = @elapsed T1, T2, ecT1, ecT2 = cas_decomposition(Casdata, refwfn.ndocc, drop_occ, actocc, actvir, ints["FOV"], ints["OOVV"], ints["OVVV"], ints["OOOV"])
    @output "Finished in {:5.5} seconds.\n" t
    Casdata = nothing

    d = [i - a for i = diag(ints["FOO"]), a = diag(ints["FVV"])]
    D = [i + j - a - b for i = diag(ints["FOO"]), j = diag(ints["FOO"]), a = diag(ints["FVV"]), b = diag(ints["FVV"])]

    RCCSD{T}(refwfn, ints, T1, T2, Alg, ecT1=ecT1, ecT2=ecT2)
end

function get_cas_data(cas::Fermi.ConfigurationInteraction.CASCI)
    
    # This function process the CAS wave function to return chunks of CI coefficients and determinants
    # corresponding to each excitation level (1, 2, 3 or 4)
    # CI coefficients are intermediate normalized.
    # The reference determinant is taken as the HF one.

    dets = cas.dets
    Ccas = cas.coef

    ref = Determinant(repeat("1",cas.ref.ndocc),repeat("1",cas.ref.ndocc))
    z = 0
    acf = 0

    @output "   • CAS Composition\n"
    @output "Excitation      N of dets\n"
    while acf < length(dets)
        x = count(d->excitation_level(ref,d)==z, dets)
        @output "{}           {}\n" z x
        z += 1
        acf += x
    end

    #Fermi.ConfigurationInteraction.showdet(ref)

    # The reference determinant is takan as the dominant configuration
    # If the reference is not the RHF determinant the C matrix has to be modified
    # such that all orbitals occupied in the new reference comes first.
    ref = dets[1]
    zeroth = repeat('1', cas.ref.ndocc)*repeat('0', cas.ref.nvir)
    hf = Fermi.ConfigurationInteraction.Determinant(zeroth, zeroth)
    if ref == hf
        @output "Dominant configuration is the RHF determinant.\n"
    else
        if ref.α ≢ ref.β
            error("Dominant determinant is not spin symmetric.")
        end
        @output "!! Dominant configuration is NOT the RHF determinant. Adapting hole-particle spaces\n"
        holes = Fermi.ConfigurationInteraction.αexclusive(hf, ref)
        particles = Fermi.ConfigurationInteraction.αexclusive(ref, hf)
        for (h,p) in zip(holes,particles)
            @inbounds for k = 1:size(cas.ref.C,1)
               cas.ref.C[k,h], cas.ref.C[k,p] = cas.ref.C[k,p], cas.ref.C[k,h]
            end

            for id in eachindex(dets)
                α = dets[id].α
                β = dets[id].β
                hbit = 1 << (h-1)
                pbit = 1 << (p-1)
                # If the bits are different, invert them
                if (α & hbit == 0) && (α & pbit != 0) 
                    α = (α | hbit) ⊻ pbit
                elseif (α & hbit != 0) && (α & pbit == 0) 
                    α = (α | pbit) ⊻ hbit
                end

                if (β & hbit == 0) && (β & pbit != 0) 
                    β = (β | hbit) ⊻ pbit
                elseif (β & hbit != 0) && (β & pbit == 0)
                    β = (β | pbit) ⊻ hbit
                end

                dets[id] = Fermi.ConfigurationInteraction.Determinant(α, β)
            end
        end
        ref = dets[1]
    end
    C0 = Ccas[1]

    # Intermediate Normalization
    abs(C0) > 1e-8 ? nothing : error("Reference coefficient is too small ($(C0)) to perform intermediate normalization")
    Ccas = Ccas ./ C0

    # Split the Cas data into excitation level
    Ccas_ex1or2 = Float64[]
    dets_ex1or2 = Determinant[]

    Ccas_ex3 = Float64[]
    dets_ex3 = Determinant[]

    Ccas_ex4 = Float64[]
    dets_ex4 = Determinant[]

    for i in eachindex(dets)

        exc = excitation_level(ref, dets[i])

        if exc == 1 || exc == 2
            push!(Ccas_ex1or2, Ccas[i])
            push!(dets_ex1or2, dets[i])

        elseif  exc == 3
            push!(Ccas_ex3, Ccas[i])
            push!(dets_ex3, dets[i])

        elseif  exc == 4
            push!(Ccas_ex4, Ccas[i])
            push!(dets_ex4, dets[i])
        end
    end

    return ref, Ccas_ex1or2, dets_ex1or2, Ccas_ex3, dets_ex3, Ccas_ex4, dets_ex4
end

#function get_casT1!(T1::Array{Float64,2}, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int)
#
#    for id in eachindex(dets)
#
#        @inbounds D = dets[id]
#        αexc = αexcitation_level(ref, D)
#        βexc = βexcitation_level(ref, D)
#
#        if αexc == 1 & βexc == 0
#            i, = αexclusive(ref, D)     
#            a, = αexclusive(D, ref) 
#
#            p = phase(ref, D)
#            # i is absolute. Take out the frozen orbitals to match the T arrays.
#            # a is abolute. Take out the occupied orbitals to get relative index.
#            @inbounds T1[i-frozen,a-ndocc] = Ccas[id]*p
#        end
#    end
#end
#
#function get_casT2!(T1::Array{Float64,2}, T2::Array{Float64,4}, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int)
#    
#    for id in eachindex(dets)
#
#        @inbounds D = dets[id]
#        αexc = αexcitation_level(ref, D)
#        βexc = βexcitation_level(ref, D)
#
#        if αexc == 1 & βexc == 1
#            elseif βexc == 1
#
#                i, = αexclusive(ref, D) 
#                j, = βexclusive(ref, D)
#                a, = αexclusive(D, ref) 
#                b, = βexclusive(D, ref) 
#
#                p = phase(ref, D)
#
#                # i and j are absolute. Take out the frozen orbitals to match the T arrays.
#                # a and b are abolute. Take out the occupied orbitals to get relative index.
#                @inbounds T2[i-frozen,j-frozen,a-ndocc,b-ndocc] = Ccas[id]*p
#            end
#
#        elseif (αexc + βexc) > 2
#            # This line relies on the fact the dets are ordered by excitation level
#            break
#        end
#    end
#
#    @tensor T2[i,j,a,b] -= T1[i,a]*T1[j,b] 
#end

function get_casT1_casT2!(T1::Array{Float64,2}, T2::Array{Float64,4}, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int)
    
    for id in eachindex(dets)

        @inbounds D = dets[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 1

            if βexc == 0

                i, = αexclusive(ref, D)     
                a, = αexclusive(D, ref) 

                p = phase(ref, D)
                # i is absolute. Take out the frozen orbitals to match the T arrays.
                # a is abolute. Take out the occupied orbitals to get relative index.
                @inbounds T1[i-frozen,a-ndocc] = Ccas[id]*p

            elseif βexc == 1

                i, = αexclusive(ref, D) 
                j, = βexclusive(ref, D)
                a, = αexclusive(D, ref) 
                b, = βexclusive(D, ref) 

                p = phase(ref, D)

                # i and j are absolute. Take out the frozen orbitals to match the T arrays.
                # a and b are abolute. Take out the occupied orbitals to get relative index.
                @inbounds T2[i-frozen,j-frozen,a-ndocc,b-ndocc] = Ccas[id]*p
            end

        elseif (αexc + βexc) > 2
            # This line relies on the fact the dets are ordered by excitation level
            break
        end
    end

    #preT2 = copy(T2)
    @tensor T2[i,j,a,b] -= T1[i,a]*T1[j,b] 
    #@tensor T12[i,j,a,b] := T1[i,a]*T1[j,b] 
    #T2 -= T12.*(T2 .!= 0)
    #for i = eachindex(T2)
    #    if T2[i] != 0.0
    #        T2[i] -= T12[i]
    #    end
    #end
    #println(preT2 ≈ T2)
end

function get_casT3!(T3::Array{Float64,4}, n::Int, f::Int, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4})

    # This function produces a particular slice of the full T3 array T3[:,:,n,:,:,f] for the αβα case

    # Clean up array
    fill!(T3, 0.0)

    for id in eachindex(dets)

        @inbounds D = dets[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 2 && βexc == 1

            # i > k, a > c
            k,i = αexclusive(ref, D)
            j,  = βexclusive(ref, D)

            if !(n in [k,i])
                continue
            end

            c,a = αexclusive(D, ref)
            b,  = βexclusive(D, ref)

            if !(f in [c,a])
                continue
            end

            # The phase is obtained by appliying annihilation and creation operations onto the referece in the order abckji
            # n and f must be in the positions k and c. a and i are the orb index that don't match n and f.

            n == k ? o1 = i : o1 = k
            f == c ? o3 = a : o3 = c

            p = 1
            _det = Determinant(ref.α, ref.β)

            _p, _det = annihilate(_det, o1, 'α')
            p = _p*p
            _p, _det = annihilate(_det, j,  'β')
            p = _p*p
            _p, _det = annihilate(_det, n,  'α')
            p = _p*p

            _p, _det = create(_det, f,  'α')
            p = _p*p
            _p, _det = create(_det, b,  'β')
            p = _p*p
            _p, _det = create(_det, o3, 'α')
            p = _p*p

            # Hole indexes are absolute. Take out the frozen orbitals to match the T arrays.
            # Particle indexes are abolute. Take out the occupied orbitals to get relative index.
            T3[o1-frozen,j-frozen,o3-ndocc,b-ndocc] = p*Ccas[id]

        end
    end

    # Shift absolute index, for cleanliness
    n = n - frozen
    f = f - ndocc

    # Arrays for decomposition
    T1_1n2f = T1[n,f]
    T1_2f = view(T1,:,f)
    T1_1n = view(T1,n,:)
    T2_1n3f = view(T2,n,:,f,:)
    T2_3f = view(T2,:,:,f,:)
    T2_1n = view(T2,n,:,:,:)
    T2_1n4f = view(T2,n,:,:,f)
    T2_2n4f = view(T2,:,n,:,f)

    @tensoropt begin
        # Decomposition of T3
        xT3[i,j,a,b] := -T1[i,a]*T1[j,b]*T1_1n2f
        xT3[i,j,a,b] +=  T1_2f[i]*T1[j,b]*T1_1n[a]
        xT3[i,j,a,b] -=  T1[i,a]*T2_1n3f[j,b]
        xT3[i,j,a,b] +=  T1_1n[a]*T2_3f[i,j,b]
        xT3[i,j,a,b] +=  T1_2f[i]*T2_1n[j,a,b] - T1_1n2f*T2[i,j,a,b]
        xT3[i,j,a,b] +=  T1[j,b]*T2_1n4f[i,a]
        xT3[i,j,a,b] -=  T1[j,b]*T2_2n4f[i,a]
    end
    for i = eachindex(T3)
        if T3[i] != 0.0
            T3[i] += xT3[i]
        end
    end
end

function get_casT4αβ!(T4::Array{Float64,4}, m::Int, n::Int, e::Int, f::Int, Ccas_ex4::Array{Float64,1}, dets_ex4::Array{Determinant,1}, Ccas_ex3::Array{Float64,1}, dets_ex3::Array{Determinant,1}, 
                   ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4}, T3_3n6f::Array{Float64,4}, T3_3m6e::Array{Float64,4})
    # Clean up array
    fill!(T4, 0.0)

    for id in eachindex(dets_ex4)

        @inbounds D = dets_ex4[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)
        if αexc == 2 && βexc == 2

            # i > k, j > l, a > c, b > d
            k,i = αexclusive(ref, D) 

            if !(m in [i,k])
                continue
            end

            l,j = βexclusive(ref, D) 

            if !(n in [l,j])
                continue
            end

            c,a = αexclusive(D,ref) 

            if !(e in [c,a])
                continue
            end

            d,b = βexclusive(D, ref) 

            if !(f in [d,b])
                continue
            end

            m == i ? o1 = k : o1 = i
            n == j ? o2 = l : o2 = j
            e == a ? o3 = c : o3 = a
            f == b ? o4 = d : o4 = b

            p = 1
            _det = Determinant(ref.α, ref.β)

            _p, _det = annihilate(_det, o1, 'α')
            p = _p*p
            _p, _det = annihilate(_det, o2, 'β')
            p = _p*p
            _p, _det = annihilate(_det, m,  'α')
            p = _p*p
            _p, _det = annihilate(_det, n,  'β')
            p = _p*p

            _p, _det = create(_det, f,  'β')
            p = _p*p
            _p, _det = create(_det, e,  'α')
            p = _p*p
            _p, _det = create(_det, o4, 'β')
            p = _p*p
            _p, _det = create(_det, o3, 'α')
            p = _p*p

            T4[o1-frozen,o2-frozen,o3-ndocc,o4-ndocc] =  p*Ccas_ex4[id]
        end
    end                 

    # Shift m,n,e, and f. Just for cleanliness
    m = m - frozen
    n = n - frozen
    e = e - ndocc
    f = f - ndocc

    #validT4 = abs.(T4) .> 10^-12

    T1_1m2e = T1[m,e]
    T1_1n2f = T1[n,f]
    pT1 = T1_1m2e*T1_1n2f
    T1_1n = view(T1,n,:)
    T1_2f = view(T1,:,f)
    T1_1m = view(T1,m,:)
    T1_2e = view(T1,:,e)

    T2_1m2n3e4f = T2[m,n,e,f]
    T2_1m3e4f = view(T2,m,:,e,f)
    T2_1n4f = view(T2,n,:,:,f)
    T2_2n4f = view(T2,:,n,:,f)
    T2_1m2n3e = view(T2,m,n,e,:)
    T2_1m3e = view(T2,m,:,e,:)
    T2_2n3e4f = view(T2,:,n,e,f)
    T2_3e4f = view(T2,:,:,e,f)
    T2_2n3e = view(T2,:,n,e,:)
    T2_3e = view(T2,:,:,e,:)
    T2_1m2n4f = view(T2,m,n,:,f)
    T2_1m4e = view(T2,m,:,:,e)
    T2_2m4e = view(T2,:,m,:,e)
    T2_1m4f = view(T2,m,:,:,f)
    T2_4f = view(T2,:,:,:,f)
    T2_1m2n = view(T2,m,n,:,:)
    T2_1m = view(T2,m,:,:,:)
    T2_2n = view(T2,:,n,:,:)

    T3_2m3n5e6f = T3_3n6f[:,m,:,e]
    T3_3n5e6f = T3_3n6f[:,:,:,e]
    T3_2m3n6f = T3_3n6f[:,m,:,:]
    T3_2n3m6e = T3_3m6e[:,n,:,:]
    T3_2n3m5f6e = T3_3m6e[:,n,:,f]
    T3_3m5f6e = T3_3m6e[:,:,:,f]

    @tensoropt begin
        xT4[i,j,a,b] := T1_1n[b]*T1[i,a]*T1_1m2e*T1_2f[j]
        xT4[i,j,a,b] -= pT1*T1[j,b]*T1[i,a]
        xT4[i,j,a,b] -= T1[j,b]*T1[i,a]*T2_1m2n3e4f
        xT4[i,j,a,b] += T1_1n[b]*T1[i,a]*T2_1m3e4f[j]
        xT4[i,j,a,b] += T1_1m2e*T1[i,a]*T2_1n4f[j,b]
        xT4[i,j,a,b] -= T1_1m2e*T1[i,a]*T2_2n4f[j,b]
        xT4[i,j,a,b] += T1_2f[j]*T1[i,a]*T2_1m2n3e[b]
        xT4[i,j,a,b] -= T1_1n2f*T1[i,a]*T2_1m3e[j,b]
        xT4[i,j,a,b] -= T1[i,a]*T3_2m3n5e6f[j,b]
        xT4[i,j,a,b] += T1[j,b]*T1_1m[a]*T1_2e[i]*T1_1n2f
        xT4[i,j,a,b] += T1[j,b]*T1_1m[a]*T2_2n3e4f[i]
        xT4[i,j,a,b] -= T1_1n[b]*T1_1m[a]*T1_2e[i]*T1_2f[j]
        xT4[i,j,a,b] -= T1_1n[b]*T1_1m[a]*T2_3e4f[i,j]
        xT4[i,j,a,b] -= T1_2e[i]*T1_1m[a]*T2_1n4f[j,b]
        xT4[i,j,a,b] += T1_2e[i]*T1_1m[a]*T2_2n4f[j,b]
        xT4[i,j,a,b] -= T1_2f[j]*T1_1m[a]*T2_2n3e[i,b]
        xT4[i,j,a,b] += T1_1n2f*T1_1m[a]*T2_3e[i,j,b]
        xT4[i,j,a,b] += T1_1m[a]*T3_3n5e6f[j,i,b]
        xT4[i,j,a,b] += T1_2e[i]*T1[j,b]*T2_1m2n4f[a]
        xT4[i,j,a,b] -= T1_1m2e*T1[j,b]*T2_2n4f[i,a]
        xT4[i,j,a,b] += T1_1n2f*T1[j,b]*T2_1m4e[i,a]
        xT4[i,j,a,b] -= T1_1n2f*T1[j,b]*T2_2m4e[i,a]
        xT4[i,j,a,b] -= T1[j,b]*T3_2n3m5f6e[i,a]
        xT4[i,j,a,b] -= T1_2e[i]*T1_1n[b]*T2_1m4f[j,a]
        xT4[i,j,a,b] += T1_1m2e*T1_1n[b]*T2_4f[i,j,a]
        xT4[i,j,a,b] -= T1_2f[j]*T1_1n[b]*T2_1m4e[i,a]
        xT4[i,j,a,b] += T1_2f[j]*T1_1n[b]*T2_2m4e[i,a]
        xT4[i,j,a,b] += T1_1n[b]*T3_3m5f6e[i,j,a]
        xT4[i,j,a,b] -= T1_2f[j]*T1_2e[i]*T2_1m2n[a,b]
        xT4[i,j,a,b] += T1_1n2f*T2_1m[j,a,b]*T1_2e[i]
        xT4[i,j,a,b] += T1_2e[i]*T3_2m3n6f[j,b,a]
        xT4[i,j,a,b] += T1_2f[j]*T1_1m2e*T2_2n[i,a,b] - pT1*T2[i,j,a,b]
        xT4[i,j,a,b] -= T1_1m2e*T3_3n6f[j,i,b,a]
        xT4[i,j,a,b] += T1_2f[j]*T3_2n3m6e[i,a,b] - T1_1n2f*T3_3m6e[i,j,a,b] - T2_1m2n3e4f*T2[i,j,a,b]
        xT4[i,j,a,b] += T2_1m3e4f[j]*T2_2n[i,a,b]
        xT4[i,j,a,b] += T2_2n3e4f[i]*T2_1m[j,a,b]
        xT4[i,j,a,b] -= T2_3e4f[i,j]*T2_1m2n[a,b]
        xT4[i,j,a,b] -= T2_1n4f[j,b]*T2_1m4e[i,a]
        xT4[i,j,a,b] += T2_1n4f[j,b]*T2_2m4e[i,a]
        xT4[i,j,a,b] += T2_2n4f[j,b]*T2_1m4e[i,a]
        xT4[i,j,a,b] -= T2_2n4f[j,b]*T2_2m4e[i,a]
        xT4[i,j,a,b] += T2_1m2n3e[b]*T2_4f[i,j,a]
        xT4[i,j,a,b] -= T2_1m3e[j,b]*T2_2n4f[i,a]
        xT4[i,j,a,b] -= T2_2n3e[i,b]*T2_1m4f[j,a]
        xT4[i,j,a,b] += T2_3e[i,j,b]*T2_1m2n4f[a]
    end

    for i = eachindex(T4)
        if T4[i] != 0.0
            T4[i] += xT4[i]
        end
    end
    
end

function get_casT4αα!(T4::Array{Float64,4}, m::Int, n::Int, e::Int, f::Int, Ccas_ex4::Array{Float64,1}, dets_ex4::Array{Determinant,1}, Ccas_ex3::Array{Float64,1}, dets_ex3::Array{Determinant,1}, 
                    ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4}, T3_3n6f::Array{Float64,4}, T3_3m6f::Array{Float64,4}, T3_3n6e::Array{Float64,4}, T3_3m6e::Array{Float64,4})
    
    fill!(T4, 0.0)

    if m == n || e == f
        return T4
    end

    for id in eachindex(dets_ex4)

        @inbounds D = dets_ex4[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 3 && βexc == 1

            # i > k > l, a > c > d
            l,k,i = αexclusive(ref, D) 

            if !(m in [k,l,i]) || !(n in [k,l,i])
                continue
            end

            j, = βexclusive(ref, D) 
            d,c,a = αexclusive(D,ref) 

            if !(e in [d,c,a]) || !(f in [d,c,a])
                continue
            end

            b, = βexclusive(D, ref) 

            o1 = filter(x-> x != m && x != n, [l,k,i])[1]
            o2 = filter(x-> x != e && x != f, [d,c,a])[1]

            p = 1
            _det = Determinant(ref.α, ref.β)

            _p, _det = annihilate(_det, o1, 'α')
            p = _p*p
            _p, _det = annihilate(_det, j,  'β')
            p = _p*p
            _p, _det = annihilate(_det, m,  'α')
            p = _p*p
            _p, _det = annihilate(_det, n,  'α')
            p = _p*p

            _p, _det = create(_det, f,  'α')
            p = _p*p
            _p, _det = create(_det, e,  'α')
            p = _p*p
            _p, _det = create(_det, b, 'β')
            p = _p*p
            _p, _det = create(_det, o2, 'α')
            p = _p*p

            T4[o1-frozen,j-frozen,o2-ndocc,b-ndocc] = p*Ccas_ex4[id]
        end
    end

    # Shift m,n,e, and f. Just for cleanliness
    m = m - frozen
    n = n - frozen
    e = e - ndocc
    f = f - ndocc
    #validT4 = T4 .!= 0.0

    # Decomposition
    T1_1m2e = T1[m,e]
    T1_1n2f = T1[n,f]
    T1_1n2e = T1[n,e]
    T1_1m2f = T1[m,f]
    T1_1m = view(T1,m,:)
    T1_2e = view(T1,:,e)
    T1_2f = view(T1,:,f)
    T1_1n = view(T1,n,:)
    pT1_a = T1_1m2e*T1_1n2f
    pT1_b = T1_1n2e*T1_1m2f

    T2_1n2m3e4f = T2[n,m,e,f]
    T2_1m2n3e4f = T2[m,n,e,f]
    T2_1n3f = view(T2,n,:,f,:)
    T2_1m3f = view(T2,m,:,f,:)
    T2_1n3e = view(T2,n,:,e,:)
    T2_1m3e = view(T2,m,:,e,:)
    T2_1n3e4f = view(T2,n,:,e,f)
    T2_2n3e4f = view(T2,:,n,e,f)
    T2_3f = view(T2,:,:,f,:)
    T2_3e = view(T2,:,:,e,:)
    T2_1m3e4f = view(T2,m,:,e,f)
    T2_2m3e4f = view(T2,:,m,e,f)
    T2_1n2m4f = view(T2,n,m,:,f)
    T2_1m2n4f = view(T2,m,n,:,f)
    T2_1n4f = view(T2,n,:,:,f)
    T2_2n4f = view(T2,:,n,:,f)
    T2_1m4f = view(T2,m,:,:,f)
    T2_2m4f = view(T2,:,m,:,f)
    T2_1n2m4e = view(T2,n,m,:,e)
    T2_1m2n4e = view(T2,m,n,:,e)
    T2_1n4e = view(T2,n,:,:,e)
    T2_2n4e = view(T2,:,n,:,e)
    T2_1m4e = view(T2,m,:,:,e)
    T2_2m4e = view(T2,:,m,:,e)
    T2_1n = view(T2,n,:,:,:)
    T2_1m = view(T2,m,:,:,:)

    T3_1m3n4e6f = view(T3_3n6f,m,:,e,:)
    T3_3n4e6f = view(T3_3n6f,:,:,e,:)
    T3_1m3n5e6f = view(T3_3n6f,m,:,:,e)
    T3_2m3n5e6f = view(T3_3n6f,:,m,:,e)
    T3_1m3n6f = view(T3_3n6f,m,:,:,:)

    T3_3m4e6f = view(T3_3m6f,:,:,e,:)
    T3_2n3m5e6f = view(T3_3m6f,:,n,:,e)

    T3_1m3n6e = view(T3_3n6e,m,:,:,:)

    @tensoropt begin
        xT4[i,j,a,b] := T1[j,b]*T1[i,a]*pT1_b
        xT4[i,j,a,b] -= T1[j,b]*T1[i,a]*pT1_a
        xT4[i,j,a,b] += T1[j,b]*T1[i,a]*T2_1n2m3e4f
        xT4[i,j,a,b] -= T1[j,b]*T1[i,a]*T2_1m2n3e4f
        xT4[i,j,a,b] -= T1_1m2e*T1[i,a]*T2_1n3f[j,b]
        xT4[i,j,a,b] += T1_1n2e*T1[i,a]*T2_1m3f[j,b]
        xT4[i,j,a,b] += T1_1m2f*T1[i,a]*T2_1n3e[j,b]
        xT4[i,j,a,b] -= T1_1n2f*T1[i,a]*T2_1m3e[j,b]
        xT4[i,j,a,b] -= T1[i,a]*T3_1m3n4e6f[j,b]
        xT4[i,j,a,b] += T1[j,b]*T1_1m[a]*T1_2e[i]*T1_1n2f
        xT4[i,j,a,b] -= T1[j,b]*T1_1m[a]*T1_1n2e*T1_2f[i]
        xT4[i,j,a,b] -= T1[j,b]*T1_1m[a]*T2_1n3e4f[i]
        xT4[i,j,a,b] += T1[j,b]*T1_1m[a]*T2_2n3e4f[i]
        xT4[i,j,a,b] += T1_2e[i]*T1_1m[a]*T2_1n3f[j,b]
        xT4[i,j,a,b] -= T1_1n2e*T1_1m[a]*T2_3f[i,j,b]
        xT4[i,j,a,b] -= T1_2f[i]*T1_1m[a]*T2_1n3e[j,b]
        xT4[i,j,a,b] += T1_1n2f*T1_1m[a]*T2_3e[i,j,b]
        xT4[i,j,a,b] += T1_1m[a]*T3_3n4e6f[i,j,b]
        xT4[i,j,a,b] -= T1[j,b]*T1_1n[a]*T1_2e[i]*T1_1m2f
        xT4[i,j,a,b] += T1[j,b]*T1_1n[a]*T1_1m2e*T1_2f[i]
        xT4[i,j,a,b] += T1[j,b]*T1_1n[a]*T2_1m3e4f[i]
        xT4[i,j,a,b] -= T1[j,b]*T1_1n[a]*T2_2m3e4f[i]
        xT4[i,j,a,b] -= T1_2e[i]*T1_1n[a]*T2_1m3f[j,b]
        xT4[i,j,a,b] += T1_1m2e*T1_1n[a]*T2_3f[i,j,b]
        xT4[i,j,a,b] += T1_2f[i]*T1_1n[a]*T2_1m3e[j,b]
        xT4[i,j,a,b] -= T1_1m2f*T1_1n[a]*T2_3e[i,j,b]
        xT4[i,j,a,b] -= T1_1n[a]*T3_3m4e6f[i,j,b]   
        xT4[i,j,a,b] -= T1_2e[i]*T1[j,b]*T2_1n2m4f[a]
        xT4[i,j,a,b] += T1_2e[i]*T1[j,b]*T2_1m2n4f[a]
        xT4[i,j,a,b] += T1_1m2e*T1[j,b]*T2_1n4f[i,a]
        xT4[i,j,a,b] -= T1_1m2e*T1[j,b]*T2_2n4f[i,a]
        xT4[i,j,a,b] -= T1_1n2e*T1[j,b]*T2_1m4f[i,a]
        xT4[i,j,a,b] += T1_1n2e*T1[j,b]*T2_2m4f[i,a]
        xT4[i,j,a,b] += T1_2f[i]*T1[j,b]*T2_1n2m4e[a]
        xT4[i,j,a,b] -= T1_2f[i]*T1[j,b]*T2_1m2n4e[a]
        xT4[i,j,a,b] -= T1_1m2f*T1[j,b]*T2_1n4e[i,a]
        xT4[i,j,a,b] += T1_1m2f*T1[j,b]*T2_2n4e[i,a]
        xT4[i,j,a,b] += T1_1n2f*T1[j,b]*T2_1m4e[i,a]
        xT4[i,j,a,b] -= T1_1n2f*T1[j,b]*T2_2m4e[i,a]
        xT4[i,j,a,b] += T1[j,b]*T3_1m3n5e6f[i,a]
        xT4[i,j,a,b] -= T1[j,b]*T3_2m3n5e6f[i,a]
        xT4[i,j,a,b] += T1[j,b]*T3_2n3m5e6f[i,a]  
        xT4[i,j,a,b] -= T1_1m2f*T1_2e[i]*T2_1n[j,a,b]
        xT4[i,j,a,b] += T1_1n2f*T1_2e[i]*T2_1m[j,a,b]
        xT4[i,j,a,b] += T1_2e[i]*T3_1m3n6f[j,a,b]
        xT4[i,j,a,b] += T1_2f[i]*T1_1m2e*T2_1n[j,a,b]
        xT4[i,j,a,b] -= pT1_a*T2[i,j,a,b]
        xT4[i,j,a,b] -= T1_1m2e*T3_3n6f[i,j,a,b]
        xT4[i,j,a,b] -= T1_2f[i]*T1_1n2e*T2_1m[j,a,b]
        xT4[i,j,a,b] += pT1_b*T2[i,j,a,b]
        xT4[i,j,a,b] += T1_1n2e*T3_3m6f[i,j,a,b]  
        xT4[i,j,a,b] -= T1_2f[i]*T3_1m3n6e[j,a,b]
        xT4[i,j,a,b] += T1_1m2f*T3_3n6e[i,j,a,b]
        xT4[i,j,a,b] -= T1_1n2f*T3_3m6e[i,j,a,b]
        xT4[i,j,a,b] += T2_1n2m3e4f*T2[i,j,a,b]
        xT4[i,j,a,b] -= T2_1m2n3e4f*T2[i,j,a,b]
        xT4[i,j,a,b] -= T2_1n3e4f[i]*T2_1m[j,a,b]
        xT4[i,j,a,b] += T2_2n3e4f[i]*T2_1m[j,a,b]
        xT4[i,j,a,b] += T2_1m3e4f[i]*T2_1n[j,a,b]
        xT4[i,j,a,b] -= T2_2m3e4f[i]*T2_1n[j,a,b]
        xT4[i,j,a,b] += T2_1n3f[j,b]*T2_1m4e[i,a]
        xT4[i,j,a,b] -= T2_1n3f[j,b]*T2_2m4e[i,a]
        xT4[i,j,a,b] -= T2_1m3f[j,b]*T2_1n4e[i,a]
        xT4[i,j,a,b] += T2_1m3f[j,b]*T2_2n4e[i,a]
        xT4[i,j,a,b] += T2_3f[i,j,b]*T2_1n2m4e[a]
        xT4[i,j,a,b] -= T2_3f[i,j,b]*T2_1m2n4e[a]
        xT4[i,j,a,b] -= T2_1n3e[j,b]*T2_1m4f[i,a]
        xT4[i,j,a,b] += T2_1n3e[j,b]*T2_2m4f[i,a]
        xT4[i,j,a,b] += T2_1m3e[j,b]*T2_1n4f[i,a]
        xT4[i,j,a,b] -= T2_1m3e[j,b]*T2_2n4f[i,a]
        xT4[i,j,a,b] -= T2_3e[i,j,b]*T2_1n2m4f[a]
        xT4[i,j,a,b] += T2_3e[i,j,b]*T2_1m2n4f[a]
    end
    for i = eachindex(T4)
        if T4[i] != 0.0
            T4[i] += xT4[i]
        end
    end
end

function get_ec_from_T3!(n::Int, f::Int, frozen::Int, ndocc::Int, ecT1::Array{Float64,2}, ecT2::Array{Float64,4}, T1::Array{Float64,2}, T3::Array{Float64,4}, fov::Array{Float64, 2}, Voovv::Array{Float64, 4}, Vovvv::Array{Float64, 4}, Vooov::Array{Float64, 4})

    # Shift m,n,e, and f. Just for cleanliness
    n = n - frozen
    f = f - ndocc

    # Arrays for ecT1 and ecT2
    Voovv_1n4f = view(Voovv, n, :, :, f)
    Voovv_2n4f = view(Voovv, :, n, :, f)
    Voovv_1n3f = view(Voovv, n, :, f, :)
    Vovvv_1n4f = view(Vovvv, n, :, :, f)
    Vovvv_1n3f = view(Vovvv, n, :, f, :)
    Vooov_1n4f = view(Vooov, n, :, :, f)
    Vooov_2n4f = view(Vooov, :, n, :, f)
    fov_1n2f = fov[n,f]
    
    @tensoropt begin
    
        # Compute ecT1
        ecT1[i,a] += 0.25*T3[m,i,e,a]*Voovv_2n4f[m,e]
        ecT1[i,a] += 1.5*T3[i,m,a,e]*Voovv_2n4f[m,e]
        ecT1[i,a] += -0.25*T3[m,i,a,e]*Voovv_2n4f[m,e]
        ecT1[i,a] += -0.5*T3[i,m,a,e]*Voovv_1n4f[m,e]
        ecT1[i,a] += -0.25*T3[m,i,e,a]*Voovv_1n4f[m,e]
        ecT1[i,a] += 0.25*T3[m,i,a,e]*Voovv_1n4f[m,e]
    
        # Compute ecT2
        ecT2[i,j,a,b] += fov_1n2f*T3[j,i,b,a]
        ecT2[i,j,a,b] += fov_1n2f*T3[i,j,a,b]
        ecT2[i,j,a,b] += -0.5*T3[i,j,e,b]*Vovvv_1n4f[a,e]
        ecT2[i,j,a,b] += 0.5*T3[i,j,e,b]*Vovvv_1n3f[a,e]
        ecT2[i,j,a,b] += T3[j,i,b,e]*Vovvv_1n3f[a,e]
        ecT2[i,j,a,b] += 0.5*T3[m,j,a,b]*Vooov_1n4f[m,i]
        ecT2[i,j,a,b] += -0.5*T3[m,j,a,b]*Vooov_2n4f[m,i]
        ecT2[i,j,a,b] -= T3[j,m,b,a]*Vooov_2n4f[m,i]
        ecT2[i,j,a,b] += 0.5*T3[m,i,b,a]*Vooov_1n4f[m,j]
        ecT2[i,j,a,b] -= T3[i,m,a,b]*Vooov_2n4f[m,j]
        ecT2[i,j,a,b] += -0.5*T3[m,i,b,a]*Vooov_2n4f[m,j]
        ecT2[i,j,a,b] += -0.5*T3[j,i,e,a]*Vovvv_1n4f[b,e]
        ecT2[i,j,a,b] += T3[i,j,a,e]*Vovvv_1n3f[b,e]
        ecT2[i,j,a,b] += 0.5*T3[j,i,e,a]*Vovvv_1n3f[b,e]
        ecT2[i,j,a,b] -= T1[m,b]*T3[i,j,a,e]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] += -0.5*T1[m,b]*T3[j,i,e,a]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] += 0.5*T1[m,b]*T3[j,i,e,a]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += 0.5*T1[m,a]*T3[i,j,e,b]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += -0.5*T1[m,a]*T3[i,j,e,b]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] -= T1[m,a]*T3[j,i,b,e]*Voovv_1n3f[m,e]
        ecT2[i,j,a,b] -= T1[j,e]*T3[i,m,a,b]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] += 0.5*T1[j,e]*T3[m,i,b,a]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += -0.5*T1[j,e]*T3[m,i,b,a]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] += 0.5*T1[i,e]*T3[m,j,a,b]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += -0.5*T1[i,e]*T3[m,j,a,b]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] -= T1[i,e]*T3[j,m,b,a]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] -= T1[m,e]*T3[i,j,a,b]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += 2.0*T1[m,e]*T3[i,j,a,b]*Voovv_2n4f[m,e]
        ecT2[i,j,a,b] -= T1[m,e]*T3[j,i,b,a]*Voovv_1n4f[m,e]
        ecT2[i,j,a,b] += 2.0*T1[m,e]*T3[j,i,b,a]*Voovv_2n4f[m,e]
    end

end

function cas_decomposition(Cas_data::Tuple, ndocc::Int, frozen::Int, actocc::Array{Int64,1}, actvir::Array{Int64,1},
                           fov::Array{Float64,2}, Voovv::Array{Float64,4}, Vovvv::Array{Float64,4}, Vooov::Array{Float64,4})

    ref, Ccas_ex1or2, dets_ex1or2, Ccas_ex3, dets_ex3, Ccas_ex4, dets_ex4 = Cas_data

    o = 1:length(actocc)
    v = 1:length(actvir)

    # Get T1 and T2
    T1 = zeros(size(fov))
    T2 = zeros(size(Voovv))
    @output "Getting T1 and T2..."
    get_casT1_casT2!(T1, T2, Ccas_ex1or2, dets_ex1or2, ref, frozen, ndocc)
    @output "Done.\n"

    # Initialize arrays
    ecT1 = zeros(size(fov))
    ecT2 = zeros(size(Voovv))

    # Allocate arrays
    T3_3n6f = similar(T2)
    T3_3m6f = similar(T2)
    T3_3n6e = similar(T2)
    T3_3m6e = similar(T2)
    #T4αβ = similar(T2)
    #T4αα = similar(T2)

    T4αβ = zeros(length(o), length(o), length(v), length(v))
    T4αα = zeros(length(o), length(o), length(v), length(v))

    # Compute ecT1
    @output "Computing ecT1, ecT2\n"
    maxn = maximum(actocc)
    for n in actocc 
        @output "{} out of {}\n" n maxn
        rn = n - frozen

        for f in actvir
            rf = f - ndocc

            get_casT3!(T3_3n6f, n, f, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)

            get_ec_from_T3!(n, f, frozen, ndocc, ecT1, ecT2, T1, T3_3n6f, fov, Voovv, Vovvv, Vooov)

            for m in actocc 
                rm = m - frozen

                get_casT3!(T3_3m6f, m, f, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)

                for e in actvir

                    get_casT3!(T3_3m6e, m, e, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)
                    get_casT3!(T3_3n6e, n, e, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)

                    get_casT4αβ!(T4αβ, m,n,e,f, Ccas_ex4, dets_ex4, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1[o,v], T2[o,o,v,v], T3_3n6f[o,o,v,v], T3_3m6e[o,o,v,v])
                    get_casT4αα!(T4αα, m,n,e,f, Ccas_ex4, dets_ex4, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1[o,v], T2[o,o,v,v], T3_3n6f[o,o,v,v], T3_3m6f[o,o,v,v], T3_3n6e[o,o,v,v], T3_3m6e[o,o,v,v])

                    re = e - ndocc
                    ecT2[o,o,v,v] += T4αβ.*Voovv[rm,rn,re,rf]
                    ecT2[o,o,v,v] += 0.25.*(T4αα + permutedims(T4αα, [2,1,4,3])).*(Voovv[rm,rn,re,rf] - Voovv[rn,rm,re,rf])
                end
            end
        end
    end

    return T1, T2, ecT1, ecT2
end
