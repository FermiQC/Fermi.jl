function Hd0(αindex::Array{Int64,1}, βindex::Array{Int64,1}, h::Array{T, 2}, V::Array{T, 4}) where T <: AbstractFloat
    """
    Σ [m|h|m] + 1/2 ΣΣ [mm|nn] - [mn|nm] 
    """
    # One electron energy
    E1 = 0.0
    # Two electron energy
    E2 = 0.0
    la = length(αindex)
    lb = length(βindex)
    for n in 1:la
        N = αindex[n]
        E1 += h[N,N]
        for m in n:la
            M = αindex[m]
            E2 += V[M,M,N,N] - V[M,N,N,M]
        end
    end
    for m in 1:la
        M = αindex[m]
        for n in m:la
            N = αindex[n]
            E2 += V[N,N,M,M] - V[N,M,M,N]
        end
        for n in 1:lb
            N = βindex[n]
            E2 += 2V[M,M,N,N]
        end
    end

    for m in 1:lb
        M = βindex[m]
        E1 += h[M,M]
        for n in m:lb
            N = βindex[n]
            E2 += V[M,M,N,N] - V[M,N,N,M]
            E2 += V[N,N,M,M] - V[N,M,M,N]
        end
    end
    return E1 + 0.5*E2
end

function Hd1(αindex::Array{Int64,1}, βindex::Array{Int64,1}, D1::Determinant, D2::Determinant, h::Array{T,2}, V::Array{T, 4}, αexc::Float64) where T <: AbstractFloat
    """
    differ m -> p
    [m|h|p] + Σ[mp|nn] - [mn|np]
    """

    # if m and p are α
    if αexc == 1
        m = first_αexclusive(D1, D2)
        p = first_αexclusive(D2, D1)

        # Compute phase by counting the number of occupied orbitals between m and p

        i = 1 << min(m,p)
        f = 1 << (max(m,p) - 1)

        ph = 1

        while i < f
            if i & D1.α ≠ 0
                ph = -ph
            end
            i = i << 1
        end

        # Compute matrix element
        E = h[m,p] 
        for i in 1:length(αindex)
            I = αindex[i]
            E += V[m,p,I,I] - V[m,I,I,p]
        end
        for i in 1:length(βindex)
            I = βindex[i]
            E += V[m,p,I,I] 
        end
        return ph*E

    else
        m = first_βexclusive(D1, D2)
        p = first_βexclusive(D2, D1)

        # Compute phase by counting the number of occupied orbitals between m and p

        i = 1 << (min(m,p))
        f = 1 << ((max(m,p) - 1))

        ph = 1

        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << (1) 
        end

        # Compute matrix element
        E = h[m,p] 
        for i in 1:length(αindex)
            I = αindex[i]
            E += V[m,p,I,I] 
        end
        for i in 1:length(βindex)
            I = βindex[i]
            E += V[m,p,I,I] - V[m,I,I,p]
        end
        return ph*E
        return E
    end
end
function Hd2(D1::Determinant, D2::Determinant, V::Array{T, 4}, αexc::Float64) where T <: AbstractFloat
    """
    mn -> pq
    [mp|nq] - [mq|np]
    """

    # If α excitation is one, it means m and n have different spins 
    if αexc == 1
        m = first_αexclusive(D1, D2)
        n = first_βexclusive(D1, D2)
        p = first_αexclusive(D2, D1)
        q = first_βexclusive(D2, D1)

        # For the phase factor, there is no interference between the two cases since they have difference spin
        # Move m <-> p
        i = 1 << (min(m,p))
        f = 1 << ((max(m,p) - 1))

        ph = 1

        while i < f
            if i & D1.α ≠ 0
                ph = -ph
            end
            i = i << (1) 
        end

        # Move n <-> q
        i = 1 << (min(n,q))
        f = 1 << ((max(n,q) - 1))

        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << (1) 
        end

        return ph*V[m,p,n,q]

    # If α excitation is two, it means m,n,p and q are all α.
    elseif αexc == 2
        m = first_αexclusive(D1, D2)
        n = second_αexclusive(D1, D2)
        p = first_αexclusive(D2, D1)
        q = second_αexclusive(D2, D1)
        
        #Transform D1 -> D2. First take m into p
        i = 1 << (min(m,p))
        f = 1 << ((max(m,p) - 1))
        
        ph = 1
        while i < f
            if i & D1.α ≠ 0
                ph = -ph
            end
            i = i << (1)
        end
        
        # Update bits
        newα = (D1.α ⊻ (1 << (m-1))) | (1 << (p-1))
        
        # Take n into q using the updated bit
        
        i = 1 << (min(n,q))
        f = 1 << ((max(n,q) - 1))
        
        while i < f
            if i & newα ≠ 0
                ph = -ph
            end
            i = i << (1&63)
        end

        return ph*(V[m,p,n,q] - V[m,q,n,p])

    # If α excitation is zero, it means m,n,p and q are all β.
    elseif αexc == 0
        m = first_βexclusive(D1, D2)
        n = second_βexclusive(D1, D2)
        p = first_βexclusive(D2, D1)
        q = second_βexclusive(D2, D1)

        #Transform D1 -> D2. First take m into p
        i = 1 << (min(m,p)&63)
        f = 1 << ((max(m,p) - 1)&63)
        
        ph = 1
        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << (1)
        end
        
        # Update bits
        newβ = (D1.β ⊻ (1 << ((m-1)))) | (1 << ((p-1)))
        
        # Take n into q using the updated bit
        
        i = 1 << (min(n,q))
        f = 1 << ((max(n,q) - 1))
        
        while i < f
            if i & newβ ≠ 0
                ph = -ph
            end
            i = i << (1)
        end

        return ph*(V[m,p,n,q] - V[m,q,n,p])
    end
end
