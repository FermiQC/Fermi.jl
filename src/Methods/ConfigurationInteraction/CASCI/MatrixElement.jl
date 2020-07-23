function Hd0(αindex::Array{Int64,1}, βindex::Array{Int64,1}, h::Array{T, 2}, V::Array{T, 4}) where T <: AbstractFloat
    """
    Σ [m|h|m] + 1/2 ΣΣ [mm|nn] - [mn|nm] 
    """
    # One electron energy
    E1 = 0.0
    # Two electron energy
    E2 = 0.0
    for m in αindex
        @inbounds E1 += h[m,m]
        for n in αindex
            @inbounds E2 += V[m,m,n,n] - V[m,n,n,m]
        end
    end

    for m in βindex
        @inbounds E1 += h[m,m]
        for n in βindex
            @inbounds E2 += V[m,m,n,n] - V[m,n,n,m]
        end
    end

    for m in αindex
        for n in βindex
            @inbounds E2 += 2V[m,m,n,n]
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
        @inbounds E = h[m,p] 
        for i in αindex
            @inbounds E += V[m,p,i,i] - V[m,i,i,p]
        end
        for i in βindex
            @inbounds E += V[m,p,i,i] 
        end
        return ph*E

    else
        m = first_βexclusive(D1, D2)
        p = first_βexclusive(D2, D1)

        # Compute phase by counting the number of occupied orbitals between m and p

        i = 1 << min(m,p)
        f = 1 << (max(m,p) - 1)

        ph = 1

        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << 1 
        end

        # Compute matrix element
        @inbounds E = h[m,p] 
        for i in αindex
            @inbounds E += V[m,p,i,i] 
        end
        for i in βindex
            @inbounds E += V[m,p,i,i] - V[m,i,i,p]
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
        i = 1 << min(m,p)
        f = 1 << (max(m,p) - 1)

        ph = 1

        while i < f
            if i & D1.α ≠ 0
                ph = -ph
            end
            i = i << 1 
        end

        # Move n <-> q
        i = 1 << min(n,q)
        f = 1 << (max(n,q) - 1)

        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << 1 
        end

        return @fastmath @inbounds ph*V[m,p,n,q]

    # If α excitation is two, it means m,n,p and q are all α.
    elseif αexc == 2
        m = first_αexclusive(D1, D2)
        n = second_αexclusive(D1, D2)
        p = first_αexclusive(D2, D1)
        q = second_αexclusive(D2, D1)
        
        #Transform D1 -> D2. First take m into p
        i = 1 << min(m,p)
        f = 1 << (max(m,p) - 1)
        
        ph = 1
        while i < f
            if i & D1.α ≠ 0
                ph = -ph
            end
            i = i << 1
        end
        
        # Update bits
        newα = (D1.α ⊻ (1 << (m-1))) | (1 << (p-1))
        
        # Take n into q using the updated bit
        
        i = 1 << min(n,q)
        f = 1 << (max(n,q) - 1)
        
        while i < f
            if i & newα ≠ 0
                ph = -ph
            end
            i = i << 1
        end

        return @fastmath @inbounds ph*(V[m,p,n,q] - V[m,q,n,p])

    # If α excitation is zero, it means m,n,p and q are all β.
    elseif αexc == 0
        m = first_βexclusive(D1, D2)
        n = second_βexclusive(D1, D2)
        p = first_βexclusive(D2, D1)
        q = second_βexclusive(D2, D1)

        #Transform D1 -> D2. First take m into p
        i = 1 << min(m,p)
        f = 1 << (max(m,p) - 1)
        
        ph = 1
        while i < f
            if i & D1.β ≠ 0
                ph = -ph
            end
            i = i << 1
        end
        
        # Update bits
        newβ = (D1.β ⊻ (1 << (m-1))) | (1 << (p-1))
        
        # Take n into q using the updated bit
        
        i = 1 << min(n,q)
        f = 1 << (max(n,q) - 1)
        
        while i < f
            if i & newβ ≠ 0
                ph = -ph
            end
            i = i << 1
        end

        return @fastmath @inbounds ph*(V[m,p,n,q] - V[m,q,n,p])
    end
end
