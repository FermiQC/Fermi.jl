function get_casT1!(T1::Array{Float64,2}, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int)

    for id in eachindex(dets)

        @inbounds D = dets[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 1 && βexc == 0

            i, = αexclusive(ref, D)     
            a, = αexclusive(D, ref) 

            p = phase(ref, D)

            # i is absolute. Take out the frozen orbitals to match the T arrays.
            # a is abolute. Take out the occupied orbitals to get relative index.
            @inbounds T1[i-frozen,a-ndocc] = Ccas[id]*p
        end
    end
end

function get_casT2!(T1::Array{Float64,2}, T2::Array{Float64,4}, Ccas::Array{Float64,1}, dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int)
    
    for id in eachindex(dets)

        @inbounds D = dets[id]
        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 1 && βexc == 1

            i, = αexclusive(ref, D) 
            j, = βexclusive(ref, D)
            a, = αexclusive(D, ref) 
            b, = βexclusive(D, ref) 

            p = phase(ref, D)

            # i and j are absolute. Take out the frozen orbitals to match the T arrays.
            # a and b are abolute. Take out the occupied orbitals to get relative index.
            @inbounds T2[i-frozen,j-frozen,a-ndocc,b-ndocc] = Ccas[id]*p - T1[i-frozen, a-ndocc]*T1[j-frozen,b-ndocc]
        end
    end
end

function get_casT3(idx::NTuple{6,Int}, Coef::Float64, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4}, ref::Determinant)

    i,j,k,a,b,c = idx

    # Create the corresponding determinant
    p = 1.0
    _p, _det = annihilate(ref, i, 'α')
    p = _p*p
    _p, _det = annihilate(_det, j,  'β')
    p = _p*p
    _p, _det = annihilate(_det, k,  'α')
    p = _p*p

    _p, _det = create(_det, c,  'α')
    p = _p*p
    _p, _det = create(_det, b,  'β')
    p = _p*p
    _p, _det = create(_det, a, 'α')
    p = _p*p

    #Cindex = findall(d->d==_det, dets) 

    #if length(Cindex) == 0
    #    return 0.0
    #end

    #T3 = Ccas[Cindex[1]]*p
    T3 = Coef*p

    # Shift absolute index, for cleanliness
    i, j, k = (i,j,k) .- frozen
    a, b, c = (a,b,c) .- ndocc

    T3 -=  T1[i,a]*T1[j,b]*T1[k,c]
    T3 +=  T1[i,c]*T1[j,b]*T1[k,a]
    T3 -=  T1[i,a]*T2[k,j,c,b]
    T3 +=  T1[k,a]*T2[i,j,c,b]
    T3 +=  T1[i,c]*T2[k,j,a,b] 
    T3 -=  T1[k,c]*T2[i,j,a,b]
    T3 +=  T1[j,b]*T2[k,i,a,c]
    T3 -=  T1[j,b]*T2[i,k,a,c]

    return T3
end

function get_casT4αβ(idx::NTuple{8,Int}, Ccas_ex4::Array{Float64,1}, dets_ex4::Array{Determinant,1}, Ccas_ex3::Array{Float64,1}, dets_ex3::Array{Determinant,1}, 
                   ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4})

    i,j,k,l,a,b,c,d = idx

    # Create the corresponding determinant
    p = 1.0

    _p, _det = annihilate(ref,  i, 'α')
    p = _p*p
    _p, _det = annihilate(_det, j, 'β')
    p = _p*p
    _p, _det = annihilate(_det, k,  'α')
    p = _p*p
    _p, _det = annihilate(_det, l,  'β')
    p = _p*p

    _p, _det = create(_det, f,  'β')
    p = _p*p
    _p, _det = create(_det, e,  'α')
    p = _p*p
    _p, _det = create(_det, o4, 'β')
    p = _p*p
    _p, _det = create(_det, o3, 'α')
    p = _p*p

    Cindex = findall(d->d==_det, dets) 

    if length(Cindex) == 0
        return 0.0
    end

    T4 =  p*Ccas_ex4[Cindex[1]]

    # Shift absolute index, for cleanliness
    i, j, k, l = (i,j,k,l) .- frozen
    a, b, c, d = (a,b,c,d) .- ndocc

    T4 -= T1[j,b]*T1[i,a]*T1[k,c]*T1[l,d]
    T4 += T1[l,b]*T1[i,a]*T1[k,c]*T1[j,d]
    T4 += T1[j,b]*T1[k,a]*T1[i,c]*T1[l,d]
    T4 -= T1[l,b]*T1[k,a]*T1[i,c]*T1[j,d]

    T4 -= T1[j,b]*T1[i,a]*T2[k,l,c,d]
    T4 += T1[l,b]*T1[i,a]*T2[k,j,c,d]
    T4 += T1[k,c]*T1[i,a]*T2[l,j,b,d]
    T4 -= T1[k,c]*T1[i,a]*T2[j,l,b,d]
    T4 += T1[j,d]*T1[i,a]*T2[k,l,c,b]
    T4 -= T1[l,d]*T1[i,a]*T2[k,j,c,b]
    T4 += T1[j,b]*T1[k,a]*T2[i,l,c,d]
    T4 -= T1[l,b]*T1[k,a]*T2[i,j,c,d]
    T4 -= T1[i,c]*T1[k,a]*T2[l,j,b,d]
    T4 += T1[i,c]*T1[k,a]*T2[j,l,b,d]
    T4 -= T1[j,d]*T1[k,a]*T2[i,l,c,b]
    T4 += T1[l,d]*T1[k,a]*T2[i,j,c,b]
    T4 += T1[i,c]*T1[j,b]*T2[k,l,a,d]
    T4 -= T1[k,c]*T1[j,b]*T2[i,l,a,d]
    T4 += T1[l,d]*T1[j,b]*T2[k,i,a,c]
    T4 -= T1[l,d]*T1[j,b]*T2[i,k,a,c]
    T4 -= T1[i,c]*T1[l,b]*T2[k,j,a,d]
    T4 += T1[k,c]*T1[l,b]*T2[i,j,a,d]
    T4 -= T1[j,d]*T1[l,b]*T2[k,i,a,c]
    T4 += T1[j,d]*T1[l,b]*T2[i,k,a,c]
    T4 -= T1[j,d]*T1[i,c]*T2[k,l,a,b]
    T4 += T1[l,d]*T1[i,c]*T2[k,j,a,b]
    T4 += T1[j,d]*T1[k,c]*T2[i,l,a,b]
    T4 -= T1[l,d]*T1[k,c]*T2[i,j,a,b]

    T4 -= T2[k,l,c,d]*T2[i,j,a,b]
    T4 += T2[k,j,c,d]*T2[i,l,a,b]
    T4 += T2[i,l,c,d]*T2[k,j,a,b]
    T4 -= T2[i,j,c,d]*T2[k,l,a,b]
    T4 -= T2[l,j,b,d]*T2[k,i,a,c]
    T4 += T2[l,j,b,d]*T2[i,k,a,c]
    T4 += T2[j,l,b,d]*T2[k,i,a,c]
    T4 -= T2[j,l,b,d]*T2[i,k,a,c]
    T4 += T2[k,l,c,b]*T2[i,j,a,d]
    T4 -= T2[k,j,c,b]*T2[i,l,a,d]
    T4 -= T2[i,l,c,b]*T2[k,j,a,d]
    T4 += T2[i,j,c,b]*T2[k,l,a,d]

    # Rename the get_casT3 function just for cleanliness
    T3(idx) = get_casT3(idx, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)

    T4 -= T1[i,a]*T3((j,k,l,b,c,d))
    T4 += T1[k,a]*T3((j,i,l,b,c,d))
    T4 -= T1[j,b]*T3((i,l,k,a,d,c))
    T4 += T1[l,b]*T3((i,j,k,a,d,c))
    T4 += T1[i,c]*T3((j,k,l,b,a,d))
    T4 -= T1[k,c]*T3((j,i,l,b,a,d))
    T4 += T1[j,d]*T3((i,l,k,a,b,c))
    T4 -= T1[l,d]*T3((i,j,k,a,b,c))

    return T4
    
end

function ec_T3onT1!(ridx::NTuple{6,Int}, T3::Float64, ecT1::Array{Float64,2}, Voovv::Array{Float64,4})

    i,j,k,a,b,c = ridx
    ecT1[j,b] += 0.25*T3*(Voovv[i,k,a,c] - Voovv[k,i,a,c])
    ecT1[i,a] += T3*(1.5*Voovv[j,k,b,c] - 0.5*Voovv[k,j,b,c])
    ecT1[j,a] += 0.25*T3*(Voovv[k,i,b,c] - Voovv[i,k,b,c])
     
end

function get_casT4αα(idx::NTuple{8,Int}, Ccas_ex4::Array{Float64,1}, dets_ex4::Array{Determinant,1}, Ccas_ex3::Array{Float64,1}, dets_ex3::Array{Determinant,1}, 
                    ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4})
    
    i,j,k,l,a,b,c,d = idx

    # Create the corresponding determinant
    p = 1.0

    _p, _det = annihilate(ref,  i, 'α')
    p = _p*p
    _p, _det = annihilate(_det, j,  'β')
    p = _p*p
    _p, _det = annihilate(_det, k,  'α')
    p = _p*p
    _p, _det = annihilate(_det, l,  'α')
    p = _p*p

    _p, _det = create(_det, a,  'α')
    p = _p*p
    _p, _det = create(_det, b,  'α')
    p = _p*p
    _p, _det = create(_det, c, 'β')
    p = _p*p
    _p, _det = create(_det, d, 'α')
    p = _p*p

    Cindex = findall(d->d==_det, dets) 

    if length(Cindex) == 0
        return 0.0
    end

    T4 =  p*Ccas_ex4[Cindex[1]]

    # Shift absolute index, for cleanliness
    i, j, k, l = (i,j,k,l) .- frozen
    a, b, c, d = (a,b,c,d) .- ndocc

    T4 -= T1[j,b]*T1[i,a]*T1[k,c]*T1[l,d]
    T4 += T1[j,b]*T1[i,a]*T1[l,c]*T1[k,d]
    T4 += T1[j,b]*T1[k,a]*T1[i,c]*T1[l,d]
    T4 -= T1[j,b]*T1[k,a]*T1[l,c]*T1[i,d]
    T4 -= T1[j,b]*T1[l,a]*T1[i,c]*T1[k,d]
    T4 += T1[j,b]*T1[l,a]*T1[k,c]*T1[i,d]

    T4 += T1[j,b]*T1[i,a]*T2[l,k,c,d]
    T4 -= T1[j,b]*T1[i,a]*T2[k,l,c,d]
    T4 -= T1[k,c]*T1[i,a]*T2[l,j,d,b]
    T4 += T1[l,c]*T1[i,a]*T2[k,j,d,b]
    T4 += T1[k,d]*T1[i,a]*T2[l,j,c,b]
    T4 -= T1[l,d]*T1[i,a]*T2[k,j,c,b]
    T4 -= T1[j,b]*T1[k,a]*T2[l,i,c,d]
    T4 += T1[j,b]*T1[k,a]*T2[i,l,c,d]
    T4 += T1[i,c]*T1[k,a]*T2[l,j,d,b]
    T4 -= T1[l,c]*T1[k,a]*T2[i,j,d,b]
    T4 -= T1[i,d]*T1[k,a]*T2[l,j,c,b]
    T4 += T1[l,d]*T1[k,a]*T2[i,j,c,b]
    T4 += T1[j,b]*T1[l,a]*T2[k,i,c,d]
    T4 -= T1[j,b]*T1[l,a]*T2[i,k,c,d]
    T4 -= T1[i,c]*T1[l,a]*T2[k,j,d,b]
    T4 += T1[k,c]*T1[l,a]*T2[i,j,d,b]
    T4 += T1[i,d]*T1[l,a]*T2[k,j,c,b]
    T4 -= T1[k,d]*T1[l,a]*T2[i,j,c,b]
    T4 -= T1[i,c]*T1[j,b]*T2[l,k,a,d]
    T4 += T1[i,c]*T1[j,b]*T2[k,l,a,d]
    T4 += T1[k,c]*T1[j,b]*T2[l,i,a,d]
    T4 -= T1[k,c]*T1[j,b]*T2[i,l,a,d]
    T4 -= T1[l,c]*T1[j,b]*T2[k,i,a,d]
    T4 += T1[l,c]*T1[j,b]*T2[i,k,a,d]
    T4 += T1[i,d]*T1[j,b]*T2[l,k,a,c]
    T4 -= T1[i,d]*T1[j,b]*T2[k,l,a,c]
    T4 -= T1[k,d]*T1[j,b]*T2[l,i,a,c]
    T4 += T1[k,d]*T1[j,b]*T2[i,l,a,c]
    T4 += T1[l,d]*T1[j,b]*T2[k,i,a,c]
    T4 -= T1[l,d]*T1[j,b]*T2[i,k,a,c]
    T4 -= T1[k,d]*T1[i,c]*T2[l,j,a,b]
    T4 += T1[l,d]*T1[i,c]*T2[k,j,a,b]
    T4 += T1[i,d]*T1[k,c]*T2[l,j,a,b]
    T4 -= T1[l,d]*T1[k,c]*T2[i,j,a,b]
    T4 -= T1[i,d]*T1[l,c]*T2[k,j,a,b]
    T4 += T1[k,d]*T1[l,c]*T2[i,j,a,b]

    T4 += T2[l,k,c,d]*T2[i,j,a,b]
    T4 -= T2[k,l,c,d]*T2[i,j,a,b]
    T4 -= T2[l,i,c,d]*T2[k,j,a,b]
    T4 += T2[i,l,c,d]*T2[k,j,a,b]
    T4 += T2[k,i,c,d]*T2[l,j,a,b]
    T4 -= T2[i,k,c,d]*T2[l,j,a,b]
    T4 += T2[l,j,d,b]*T2[k,i,a,c]
    T4 -= T2[l,j,d,b]*T2[i,k,a,c]
    T4 -= T2[k,j,d,b]*T2[l,i,a,c]
    T4 += T2[k,j,d,b]*T2[i,l,a,c]
    T4 += T2[i,j,d,b]*T2[l,k,a,c]
    T4 -= T2[i,j,d,b]*T2[k,l,a,c]
    T4 -= T2[l,j,c,b]*T2[k,i,a,d]
    T4 += T2[l,j,c,b]*T2[i,k,a,d]
    T4 += T2[k,j,c,b]*T2[l,i,a,d]
    T4 -= T2[k,j,c,b]*T2[i,l,a,d]
    T4 -= T2[i,j,c,b]*T2[l,k,a,d]
    T4 += T2[i,j,c,b]*T2[k,l,a,d]

    # Rename the get_casT3 function just for cleanliness
    T3(idx) = get_casT3(idx, Ccas_ex3, dets_ex3, ref, frozen, ndocc, T1, T2)

    T4 -= T1[i,a]*T3((k,j,l,c,b,d))
    T4 += T1[k,a]*T3((i,j,l,c,b,d))
    T4 -= T1[l,a]*T3((i,j,k,c,b,d))
    T4 += T1[j,b]*T3((k,i,l,a,c,d))
    T4 -= T1[j,b]*T3((i,k,l,a,c,d))
    T4 += T1[j,b]*T3((i,l,k,a,c,d))
    T4 += T1[i,c]*T3((k,j,l,a,b,d))
    T4 -= T1[k,c]*T3((i,j,l,a,b,d))
    T4 += T1[l,c]*T3((i,j,k,a,b,d))
    T4 -= T1[i,d]*T3((k,j,l,a,b,c))
    T4 += T1[k,d]*T3((i,j,l,a,b,c))
    T4 -= T1[l,d]*T3((i,j,k,a,b,c))

end

function new_cas_decomposition(Cas_data::Tuple, ndocc::Int, frozen::Int, fov::Array{Float64,2}, Voovv::Array{Float64,4}, Vovvv::Array{Float64,4}, Vooov::Array{Float64,4})

    ref, Ccas_ex1or2, dets_ex1or2, Ccas_ex3, dets_ex3, Ccas_ex4, dets_ex4 = Cas_data

    # Get T1 and T2
    T1 = zeros(size(fov))
    T2 = zeros(size(Voovv))
    @output "Getting T1..."
    get_casT1!(T1, Ccas_ex1or2, dets_ex1or2, ref, frozen, ndocc)
    @output " Done.\n"
    @output "Getting T2..."
    get_casT2!(T1, T2, Ccas_ex1or2, dets_ex1or2, ref, frozen, ndocc)
    @output " Done.\n"

    # Initialize arrays
    ecT1 = zeros(size(fov))
    ecT2 = zeros(size(Voovv))

    # Compute ecT1
    for n in eachindex(dets_ex3)
        D = dets_ex3[n]
        C = Ccas_ex3[n]

        if abs(C) < 10^-10
            continue
        end

        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc != 2 || βexc != 1
            continue
        end

        k,i = αexclusive(ref, D)
        j,  = βexclusive(ref,D)
        c,a = αexclusive(D, ref)
        b,  = βexclusive(D, ref)

        T3 = get_casT3((i,j,k,a,b,c), C, frozen, ndocc, T1, T2, ref)
        i, j, k = (i,j,k) .- frozen
        a, b, c = (a,b,c) .- ndocc

        ec_T3onT1!((i,j,k,a,b,c), T3, ecT1, Voovv)
        ec_T3onT1!((k,j,i,a,b,c),-T3, ecT1, Voovv)
        ec_T3onT1!((i,j,k,c,b,a),-T3, ecT1, Voovv)
        ec_T3onT1!((k,j,i,c,b,a), T3, ecT1, Voovv)
    end

    return ecT1
end
