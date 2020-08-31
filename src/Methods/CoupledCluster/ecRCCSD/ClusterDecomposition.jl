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

    T3 = Coef*p

    # Shift absolute index, for cleanliness
    i, j, k = (i,j,k) .- frozen
    a, b, c = (a,b,c) .- ndocc

    @inbounds begin
        T3 -=  T1[i,a]*T1[j,b]*T1[k,c]
        T3 +=  T1[i,c]*T1[j,b]*T1[k,a]
        T3 -=  T1[i,a]*T2[k,j,c,b]
        T3 +=  T1[k,a]*T2[i,j,c,b]
        T3 +=  T1[i,c]*T2[k,j,a,b] 
        T3 -=  T1[k,c]*T2[i,j,a,b]
        T3 +=  T1[j,b]*T2[k,i,a,c]
        T3 -=  T1[j,b]*T2[i,k,a,c]
    end
    return T3
end

function find_casT3(ridx::NTuple{6,Int}, C3coef::Array{Float64,1}, C3dets::Array{Determinant,1}, ref::Determinant, frozen::Int, ndocc::Int, T1::Array{Float64,2}, T2::Array{Float64,4})

    i,j,k,a,b,c = ridx
    # Change indexes to absolute
    i, j, k = (i,j,k) .+ frozen
    a, b, c = (a,b,c) .+ ndocc

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

    search = findall(d->d==_det, C3dets)

    if length(search) == 0
        return 0.0
    end

    Coef = C3coef[search[1]]
    T3 = Coef*p

    # Shift absolute index, for cleanliness
    i, j, k = (i,j,k) .- frozen
    a, b, c = (a,b,c) .- ndocc

    @inbounds begin
        T3 -=  T1[i,a]*T1[j,b]*T1[k,c]
        T3 +=  T1[i,c]*T1[j,b]*T1[k,a]
        T3 -=  T1[i,a]*T2[k,j,c,b]
        T3 +=  T1[k,a]*T2[i,j,c,b]
        T3 +=  T1[i,c]*T2[k,j,a,b] 
        T3 -=  T1[k,c]*T2[i,j,a,b]
        T3 +=  T1[j,b]*T2[k,i,a,c]
        T3 -=  T1[j,b]*T2[i,k,a,c]
    end
    return T3
end

function get_casT4αβ(idx::NTuple{8,Int}, Coef::Float64, C3coef::Array{Float64,1}, C3dets::Array{Determinant,1}, 
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

    _p, _det = create(_det, d, 'β')
    p = _p*p
    _p, _det = create(_det, c, 'α')
    p = _p*p
    _p, _det = create(_det, b, 'β')
    p = _p*p
    _p, _det = create(_det, a, 'α')
    p = _p*p

    T4 =  p*Coef

    # Shift absolute index, for cleanliness
    i, j, k, l = (i,j,k,l) .- frozen
    a, b, c, d = (a,b,c,d) .- ndocc

    @inbounds begin
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
    end

    # Rename the find_T3 function just for cleanliness
    T3(ridx) = find_casT3(ridx, C3coef, C3dets, ref, frozen, ndocc, T1, T2)

    @inbounds begin
        T4 -= T1[i,a]*T3((j,k,l,b,c,d))
        T4 += T1[k,a]*T3((j,i,l,b,c,d))
        T4 -= T1[j,b]*T3((i,l,k,a,d,c))
        T4 += T1[l,b]*T3((i,j,k,a,d,c))
        T4 += T1[i,c]*T3((j,k,l,b,a,d))
        T4 -= T1[k,c]*T3((j,i,l,b,a,d))
        T4 += T1[j,d]*T3((i,l,k,a,b,c))
        T4 -= T1[l,d]*T3((i,j,k,a,b,c))
    end
    return T4
end

function get_casT4αα(idx::NTuple{8,Int}, Coef::Float64, C3coef::Array{Float64,1}, C3dets::Array{Determinant,1}, 
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

    _p, _det = create(_det, d,  'α')
    p = _p*p
    _p, _det = create(_det, c,  'α')
    p = _p*p
    _p, _det = create(_det, b, 'β')
    p = _p*p
    _p, _det = create(_det, a, 'α')
    p = _p*p

    T4 =  p*Coef

    # Shift absolute index, for cleanliness
    i, j, k, l = (i,j,k,l) .- frozen
    a, b, c, d = (a,b,c,d) .- ndocc

    @inbounds begin
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
    end

    # Rename the find_T3 function just for cleanliness
    T3(ridx) = find_casT3(ridx, C3coef, C3dets, ref, frozen, ndocc, T1, T2)

    @inbounds begin
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
    return T4
end

function ec_T3onT1!(ridx::NTuple{6,Int}, T3::Float64, ecT1::Array{Float64,2}, Voovv::Array{Float64,4})

    i,j,k,a,b,c = ridx
    @inbounds begin
        ecT1[j,b] += 0.25*T3*(Voovv[i,k,a,c] - Voovv[k,i,a,c])
        ecT1[i,a] += T3*(1.5*Voovv[j,k,b,c] - 0.5*Voovv[k,j,b,c])
        ecT1[j,a] += 0.25*T3*(Voovv[k,i,b,c] - Voovv[i,k,b,c])
    end
end

function ec_T3onT2!(ridx::NTuple{6,Int}, T3::Float64, ecT2::Array{Float64,4}, T1::Array{Float64,2}, fov::Array{Float64,2}, Voovv::Array{Float64,4}, Vovvv::Array{Float64,4}, Vooov::Array{Float64,4})

    i,j,k,a,b,c = ridx

    # Simple
    @inbounds begin
        ecT2[j,i,b,a] += fov[k,c]*T3
        ecT2[i,j,a,b] += fov[k,c]*T3

        # Slice
        ecT2[j,i,:,a] += T3*Vovvv[k,:,c,b]
        ecT2[:,j,a,b] += 0.5*T3*Vooov[k,i,:,c]
        ecT2[:,j,a,b] += -0.5*T3*Vooov[i,k,:,c]
        ecT2[:,i,b,a] -= T3*Vooov[j,k,:,c]
        ecT2[j,:,b,a] += 0.5*T3*Vooov[k,i,:,c]
        ecT2[i,:,a,b] -= T3*Vooov[j,k,:,c]
        ecT2[j,:,b,a] += -0.5*T3*Vooov[i,k,:,c]
        ecT2[j,i,b,:] += -0.5*T3*Vovvv[k,:,a,c]
        ecT2[i,j,a,:] += T3*Vovvv[k,:,c,b]
        ecT2[j,i,b,:] += 0.5*T3*Vovvv[k,:,c,a]
        ecT2[i,j,:,b] += 0.5*T3*Vovvv[k,:,c,a]
        ecT2[i,j,:,b] += -0.5*T3*Vovvv[k,:,a,c]

        # Contraction
        sI,sA = size(T1)
        for I = 1:sI
            for A = 1:sA
                ecT2[j,i,A,a] -= T1[I,A]*T3*Voovv[k,I,c,b]
                ecT2[j,i,b,A] += 0.5*T1[I,A]*T3*Voovv[k,I,a,c]
                ecT2[j,i,b,A] += -0.5*T1[I,A]*T3*Voovv[I,k,a,c]
                ecT2[i,j,a,A] -= T1[I,A]*T3*Voovv[I,k,b,c]
                ecT2[i,I,a,b] -= T1[I,A]*T3*Voovv[j,k,A,c]
                ecT2[j,I,b,a] += 0.5*T1[I,A]*T3*Voovv[k,i,A,c]
                ecT2[j,I,b,a] += -0.5*T1[I,A]*T3*Voovv[i,k,A,c]
                ecT2[I,j,a,b] += 0.5*T1[I,A]*T3*Voovv[k,i,A,c]
                ecT2[I,j,a,b] += -0.5*T1[I,A]*T3*Voovv[i,k,A,c]
                ecT2[I,i,b,a] -= T1[I,A]*T3*Voovv[j,k,A,c]
                ecT2[i,j,a,b] -= T1[I,A]*T3*Voovv[k,I,A,c]
                ecT2[i,j,a,b] += 2.0*T1[I,A]*T3*Voovv[I,k,A,c]
                ecT2[j,i,b,a] -= T1[I,A]*T3*Voovv[k,I,A,c]
                ecT2[j,i,b,a] += 2.0*T1[I,A]*T3*Voovv[I,k,A,c]
                ecT2[i,j,A,b] += -0.5*T1[I,A]*T3*Voovv[I,k,a,c]
                ecT2[i,j,A,b] += 0.5*T1[I,A]*T3*Voovv[k,I,a,c]
            end
        end
    end
end

function ec_T4αβonT2!(ridx::NTuple{8,Int}, T4::Float64, ecT2::Array{Float64,4}, Voovv::Array{Float64,4})

    i,j,k,l,a,b,c,d = ridx

    ecT2[i,j,a,b] += T4*Voovv[l,k,d,c]
end

function ec_T4ααonT2!(ridx::NTuple{8,Int}, T4::Float64, ecT2::Array{Float64,4}, Voovv::Array{Float64,4})

    i,j,k,l,a,b,c,d = ridx

    ec = 0.25*T4*(Voovv[k,l,c,d] - Voovv[l,k,c,d])
    ecT2[i,j,a,b] += ec
    ecT2[j,i,b,a] += ec
end

function new_cas_decomposition(Cas_data::Tuple, ndocc::Int, frozen::Int, fov::Array{Float64,2}, Voovv::Array{Float64,4}, Vovvv::Array{Float64,4}, Vooov::Array{Float64,4})

    ref, Ccas_ex1or2, dets_ex1or2, C3coef, C3dets, C4coef, C4dets = Cas_data

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

    @output "Computing external correction from T3..."
    # Compute External Correction from T3
    for n in eachindex(C3dets)
        D = C3dets[n]
        C = C3coef[n]

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

        # Include all permutations of T3 into ecT1
        ec_T3onT1!((i,j,k,a,b,c), T3, ecT1, Voovv)
        ec_T3onT1!((k,j,i,a,b,c),-T3, ecT1, Voovv)
        ec_T3onT1!((i,j,k,c,b,a),-T3, ecT1, Voovv)
        ec_T3onT1!((k,j,i,c,b,a), T3, ecT1, Voovv)

        # Include all permutations of T3 into ecT2
        ec_T3onT2!((i,j,k,a,b,c), T3, ecT2, T1, fov, Voovv, Vovvv, Vooov)
        ec_T3onT2!((k,j,i,a,b,c),-T3, ecT2, T1, fov, Voovv, Vovvv, Vooov)
        ec_T3onT2!((i,j,k,c,b,a),-T3, ecT2, T1, fov, Voovv, Vovvv, Vooov)
        ec_T3onT2!((k,j,i,c,b,a), T3, ecT2, T1, fov, Voovv, Vovvv, Vooov)
    end
    @output "Done\n"

    # Compute External Correction from T4
    @output "Computing external correction from T4..."
    for n in eachindex(C4dets)
        D = C4dets[n]
        C = C4coef[n]

        if abs(C) < 10^-10
            continue
        end

        αexc = αexcitation_level(ref, D)
        βexc = βexcitation_level(ref, D)

        if αexc == 2 && βexc == 2

            k,i  = αexclusive(ref, D)
            l,j  = βexclusive(ref,D)
            c,a  = αexclusive(D, ref)
            d,b  = βexclusive(D, ref)

            T4 = get_casT4αβ((i,j,k,l,a,b,c,d), C, C3coef, C3dets, ref, frozen, ndocc, T1, T2)

            i, j, k, l = (i,j,k,l) .- frozen
            a, b, c, d = (a,b,c,d) .- ndocc

            # Include all permutations of T4αβ into ecT2
            ec_T4αβonT2!((i,j,k,l,a,b,c,d),  T4, ecT2, Voovv)
            ec_T4αβonT2!((i,j,k,l,a,d,c,b), -T4, ecT2, Voovv)
            ec_T4αβonT2!((i,j,k,l,c,d,a,b),  T4, ecT2, Voovv)
            ec_T4αβonT2!((i,j,k,l,c,b,a,d), -T4, ecT2, Voovv)
            ec_T4αβonT2!((k,j,i,l,a,d,c,b),  T4, ecT2, Voovv)
            ec_T4αβonT2!((k,j,i,l,a,b,c,d), -T4, ecT2, Voovv)
            ec_T4αβonT2!((k,j,i,l,c,d,a,b), -T4, ecT2, Voovv)
            ec_T4αβonT2!((k,j,i,l,c,b,a,d),  T4, ecT2, Voovv)
            ec_T4αβonT2!((k,l,i,j,a,d,c,b), -T4, ecT2, Voovv)
            ec_T4αβonT2!((k,l,i,j,a,b,c,d),  T4, ecT2, Voovv)
            ec_T4αβonT2!((k,l,i,j,c,d,a,b),  T4, ecT2, Voovv)
            ec_T4αβonT2!((k,l,i,j,c,b,a,d), -T4, ecT2, Voovv)
            ec_T4αβonT2!((i,l,k,j,a,d,c,b),  T4, ecT2, Voovv)
            ec_T4αβonT2!((i,l,k,j,a,b,c,d), -T4, ecT2, Voovv)
            ec_T4αβonT2!((i,l,k,j,c,d,a,b), -T4, ecT2, Voovv)
            ec_T4αβonT2!((i,l,k,j,c,b,a,d),  T4, ecT2, Voovv)

        elseif αexc == 3 && βexc == 1

            l,k,i  = αexclusive(ref, D)
            j,     = βexclusive(ref,D)
            d,c,a  = αexclusive(D, ref)
            b,     = βexclusive(D, ref)

            T4 = get_casT4αα((i,j,k,l,a,b,c,d), C, C3coef, C3dets, ref, frozen, ndocc, T1, T2)

            i, j, k, l = (i,j,k,l) .- frozen
            a, b, c, d = (a,b,c,d) .- ndocc

            # Include all permutations of T4αβ into ecT2
            ec_T4ααonT2!((i,j,k,l,a,b,c,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,k,l,c,b,d,a),  T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,k,l,a,b,d,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,k,l,c,b,a,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,k,l,d,b,a,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,k,l,d,b,c,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,c,b,d,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,a,b,d,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,c,b,a,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,a,b,c,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,d,b,a,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,i,l,d,b,c,a),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,c,b,d,a),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,a,b,d,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,c,b,a,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,a,b,c,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,d,b,a,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((k,j,l,i,d,b,c,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,c,b,d,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,a,b,d,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,c,b,a,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,a,b,c,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,d,b,a,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,k,i,d,b,c,a),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,c,b,d,a),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,a,b,d,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,c,b,a,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,a,b,c,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,d,b,a,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((l,j,i,k,d,b,c,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,c,b,d,a), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,a,b,d,c),  T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,c,b,a,d),  T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,a,b,c,d), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,d,b,a,c), -T4, ecT2, Voovv)
            ec_T4ααonT2!((i,j,l,k,d,b,c,a),  T4, ecT2, Voovv)
        end
    end
    @output "Done\n"
    return T1, T2, ecT1, ecT2
end
