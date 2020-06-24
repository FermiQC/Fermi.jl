module MatrixElement
using TensorOperations
using Fermi.ConfigurationInteraction.DetOperations
using LinearAlgebra

export Hd0
export Hd1
export Hd2

"""
Module implementing matrix elements (currently just Hamiltonian) via
Slater's rules

Formulas and notation from Szabo and Ostlund p. 70

syntax --> <nparticle>Electron_<matrix>d<ndiff>
	   -->
	   --> e.g. for the one particle matrix element of a Hamiltonian betweeen
	   --> two determinants differing by 1 orbital:
	   --> 	OneElectron_Hd1
"""

using Fermi.Wavefunction
using Fermi.ConfigurationInteraction.DetOperations

export Hd0
export Hd1
export Hd2

function Hd0(αindex::Array{Int64,1}, βindex::Array{Int64,1}, h::Array{Float64, 2}, V::Array{Float64, 4})
    """
    Σ [m|h|m] + 1/2 ΣΣ [mm|nn] - [mn|nm] 
    """

    # One-electron contribution
    E = tr(h[αindex, αindex]) + tr(h[βindex, βindex])

    # Two-electron contributions

    _Vα  = V[αindex, αindex, αindex, αindex]
    _Vβ  = V[βindex, βindex, βindex, βindex]
    _Vαβ = V[αindex, αindex, βindex, βindex]

    @tensor E += 0.5*(_Vα[m,m,n,n] + _Vβ[m,m,n,n] + 2_Vαβ[m,m,n,n] - _Vα[m,n,n,m] - _Vβ[m,n,n,m])

    return E

end

function Hd1(αindex::Array{Int64,1}, βindex::Array{Int64,1}, D1::Determinant, D2::Determinant, h::Array{Float64,2}, V::Array{Float64, 4}, αexc::Float64)
    """
    differ m -> p
    [m|h|p] + Σ[mp|nn] - [mn|np]
    """

    # if m and p are α
    if αexc == 1
        m, = αexclusive(D1, D2)
        p, = αexclusive(D2, D1)

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

        # αindex and βindex are use to mask the two electron integrals such that we can take the trace
        # over the occupied electrons. The two electron terms are, respectively: Jα, Jβ, and Kα
        @inbounds E = ph*(h[m,p] + tr(V[m, p, αindex, αindex]) + tr(V[m, p, βindex, βindex]) - tr(V[m, αindex, αindex, p]))
        return E

    else
        m, = βexclusive(D1, D2)
        p, = βexclusive(D2, D1)

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

        # αindex and βindex are use to mask the two electron integrals such that we can take the trace
        # over the occupied electrons. The two electron terms are, respectively: Jα, Jβ, and Kα
        @inbounds E = ph*(h[m,p] + tr(V[m, p, αindex, αindex]) + tr(V[m, p, βindex, βindex]) - tr(V[m, βindex, βindex, p]))
        return E
    end
end

function Hd2(D1::Determinant, D2::Determinant, V::Array{Float64, 4}, αexc::Float64)
    """
    mn -> pq
    [mp|nq] - [mq|np]
    """

    # If α excitation is one, it means m and n have different spins 
    if αexc == 1
        m, = αexclusive(D1, D2)
        n, = βexclusive(D1, D2)
        p, = αexclusive(D2, D1)
        q, = βexclusive(D2, D1)

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
        m,n = αexclusive(D1, D2)
        p,q = αexclusive(D2, D1)

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
        m,n = βexclusive(D1, D2)
        p,q = βexclusive(D2, D1)

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

end #module
