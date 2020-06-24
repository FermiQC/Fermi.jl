"""
    Fermi.ConfigurationInteraction.DetOperations

Module containing Determinant objects and associated operations.

## Structures

    Determinant -> Holds alpha and beta strings for a determinant as integers

## Functions

    αlist -> Return the alpha string of a Determinant as a list
    βlist -> Return the beta string of a Determinant as a list
    αindex -> Return the indexes of the occupied alpha electrons in a Determinant
    βindex -> Return the indexes of the occupied beta electrons in a Determinant
    αexcitation_level -> Compare two determinants to return the excitation level of the alpha electrons
    βexcitation_level -> Compare two determinants to return the excitation level of the beta electrons
    excitation_level ->  Compare two determinants to return the excitation level
    αexclusive -> Compare two determinants to return the indexes of the alpha electrons present in 
    the first, but not in the second
    βexclusive -> Compare two determinants to return the indexes of the beta electrons present in 
    the first, but not in the second
    exclusive -> Returns a list with tuples (orbital index, spin) of orbitals populated in the first, but
    not in the second
    annihilate -> Creates a copy of the determinat where the specified electrons was deleted. Returns the 
    phase factor along with the new determinant
    create -> Creates a copy of the determinat where the an electron was added to a specified orbital. 
    Returns the phase factor along with the new determinant
    phase -> Returns the phase (+1 or -1) associated with the transformation of one determinant
    into the other through second quantization operators
    showdet -> Prints alpha and beta strings of a Determinant

"""
module DetOperations

export Determinant
export αlist
export βlist
export αindex   
export βindex  
export αexcitation_level
export βexcitation_level
export excitation_level
export αexclusive
export βexclusive
export exclusive
export annihilate
export create
export phase
export showdet

"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant

    struct -> holding alpha and beta strings for a determinant as integers.

## Fields

    α::Int -> integer representing an alpha string
    β::Int -> integer representing an beta string
"""
struct Determinant
    α::Int
    β::Int
end

"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant(α::String, β::String)

    Constructor function for Determinant object using strings.

## Arguments

    α::String -> Alpha string ordered from left to right   
    β::String -> Beta string ordered from left to right
"""
function Determinant(α::String, β::String)
    
    αint = parse(Int, reverse(α); base=2) 
    βint = parse(Int, reverse(β); base=2) 

    Determinant(αint, βint)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αlist(D::Determinant)

    Return the alpha string of a Determinant as a list

## Arguments

    D::Determinant -> Determinant which the alpha list is to be extracted
"""
function αlist(D::Determinant)

    return [parse(Int, ss) for ss in reverse(bitstring(D.α))]
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βlist(D::Determinant)

    Return the beta string of a Determinant as a list

## Arguments

    D::Determinant -> Determinant which the beta list is to be extracted
"""
function βlist(D::Determinant)

    return [parse(Int, ss) for ss in reverse(bitstring(D.β))]
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αindex(D::Determinant)

    Return the indexes of the occupied alpha electrons in a Determinant

## Arguments

    D::Determinant -> Determinant which the indexes are to be extracted
    N::Int -> Number of electron indexes to be returned
"""
function αindex(D::Determinant, N::Int)

    out = Array{Int64,1}(undef,N)
    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ N
        if 1<<(i-1) & D.α ≠ 0
            out[e] = i
            e += 1
        end
        i += 1
    end

    return out

end

"""
    Fermi.ConfigurationInteraction.DetOperations.βindex(D::Determinant)

    Return the indexes of the occupied beta electrons in a Determinant

## Arguments

    D::Determinant -> Determinant which the indexes are to be extracted
    N::Int -> Number of electron indexes to be returned
"""
function βindex(D::Determinant, N::Int)

    out = Array{Int64,1}(undef,N)
    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ N
        if 1<<(i-1) & D.β ≠ 0
            out[e] = i
            e += 1
        end
        i += 1
    end

    return out

end

"""
    Fermi.ConfigurationInteraction.DetOperations.αexcitation_level(D1::Determinant, D2::Determinant)

    Compare two determinants to return the excitation level of the alpha electrons

## Arguments

    D1::Determinant -> Determinant to be compared
    D2::Determinant -> Determinant to be compared
"""
function αexcitation_level(D1::Determinant, D2::Determinant)

    αdiff = D1.α ⊻ D2.α
    exc = 0
    i = 1
    while i ≤ αdiff
        if i & αdiff ≠ 0
            exc += 1
        end
        i = i << 1 
    end

    return exc/2
        
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexcitation_level(D1::Determinant, D2::Determinant)

    Compare two determinants to return the excitation level of the beta electrons

## Arguments

    D1::Determinant -> Determinant to be compared
    D2::Determinant -> Determinant to be compared
"""
function βexcitation_level(D1::Determinant, D2::Determinant)

    βdiff = D1.β ⊻ D2.β
    exc = 0
    i = 1
    while i ≤ βdiff
        if i & βdiff ≠ 0
            exc += 1
        end
        i = i << 1 
    end

    return exc/2

end

"""
    Fermi.ConfigurationInteraction.DetOperations.excitation_level(D1::Determinant, D2::Determinant)

    Compare two determinants to return the excitation level

## Arguments

    D1::Determinant -> Determinant to be compared
    D2::Determinant -> Determinant to be compared
"""
function excitation_level(D1::Determinant, D2::Determinant)

    return αexcitation_level(D1, D2) + βexcitation_level(D1,D2)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αexclusive(D1::Determinant, D2::Determinant)

    Compare two determinants to return the indexes of the alpha electrons present in the first, but
    not in the second

## Arguments

    D1::Determinant -> Determinant where the alpha electrons must be
    D2::Determinant -> Determinant where the alpha electrons must not be
"""
function αexclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while 1<<(i-1) ≤ αexcl
        if 1<<(i-1) & αexcl ≠ 0
            push!(out, i)
        end
        i += 1
    end

    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexclusive(D1::Determinant, D2::Determinant)

    Compare two determinants to return the indexes of the beta electrons present in the first, but
    not in the second

## Arguments

    D1::Determinant -> Determinant where the beta electrons must be
    D2::Determinant -> Determinant where the beta electrons must not be
"""
function βexclusive(D1::Determinant, D2::Determinant)

    βexcl = D1.β ⊻ D2.β & D1.β

    out = []
    i = 1
    # Save betas exclusives, in crescent order
    while 1<<(i-1) ≤ βexcl
        if 1<<(i-1) & βexcl ≠ 0
            push!(out, i)
        end
        i += 1
    end

    return out

end

"""
    Fermi.ConfigurationInteraction.DetOperations.exclusive(D1::Determinant, D2::Determinant)

    Returns a list with tuples (orbital index, spin) of orbitals populated in the first, but
    not in the second

## Arguments

    D1::Determinant -> Determinant where the electrons must be
    D2::Determinant -> Determinant where the electrons must not be
"""
function exclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α
    βexcl = D1.β ⊻ D2.β & D1.β

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while 1<<(i-1) ≤ αexcl
        if 1<<(i-1) & αexcl ≠ 0
            push!(out, (i, 'α'))
        end
        i += 1
    end
    i = 1
    # Save betas exclusives, in crescent order
    while 1<<(i-1) ≤ βexcl
        if 1<<(i-1) & βexcl ≠ 0
            push!(out, (i, 'β'))
        end
        i += 1
    end
    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.annihilate(D::Determinant, orb::Int, spin::Char)

    Creates a copy of the determinat where the specified electrons was deleted. Returns the phase
    factor along with the new determinant

## Arguments

    D::Determinant -> Determinant where the electron will be annihilate
    orb::Int -> Index of the orbital where the electron will be annihilated
    spin::Char -> Spin of the electron, must be 'α' or 'β'
"""
function annihilate(D::Determinant, orb::Int, spin::Char)

    if spin == 'α'
        #if D.α & (1 << (orb-1)) == 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = 0
        i = 1
        while i < (1 << (orb-1))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α ⊻ (1 << (orb-1))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        #if D.β & (1 << (orb-1)) == 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = 1
        while i < (1 << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β ⊻ (1 << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

"""
    Fermi.ConfigurationInteraction.DetOperations.create(D::Determinant, orb::Int, spin::Char)

    Creates a copy of the determinat where the an electron was added to a specified orbital. 
    Returns the phase factor along with the new determinant

## Arguments

    D::Determinant -> Determinant where the electron will be created
    orb::Int -> Index of the orbital where the electron will be created
    spin::Char -> Spin of the electron, must be 'α' or 'β'
"""
function create(D::Determinant, orb::Int, spin::Char)

    if spin == 'α'
        #if D.α & (1 << (orb-1)) ≠ 0
        #    error("Creation error. Orbital $orb is occupied")
        #end

        # Determine sign
        l = 0
        i = 1
        while i < (1 << (orb-1))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α | (1 << (orb-1))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        #if D.β & (1 << (orb-1)) ≠ 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = 1
        while i < (1 << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β | (1 << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

#function α1phase

"""
    Fermi.ConfigurationInteraction.DetOperations.phase(D1::Determinant, D2::Determinant)

    Returns the phase (+1 or -1) associated with the transformation of one determinant
    into the other through second quantization operators

## Arguments

    D1::Determinant -> Starting determinant
    D2::Determinant -> Final determinant
"""
function phase(D1::Determinant, D2::Determinant)

    p = 1
    _det = Determinant(D1.α, D1.β)

    # For a string of excitation operators: abc...kji. We apply the annihilation operations such that k < j < i
    # Thus, the reverse function
    for (i,σ) in reverse(exclusive(D1, D2))
        f, _det = annihilate(_det, i, σ)
        p = f*p
    end
    # For the creation operations a > b > c. Thus, no reverse.
    for (a,σ) in exclusive(D2, D1)
        f, _det = create(_det, a, σ)
        p = f*p
    end

    return p
end

"""
    Fermi.ConfigurationInteraction.DetOperations.phase(D::Determinant, l::Int=0)

    Prints alpha and beta strings of a Determinant

## Arguments

    D::Determinant -> Determinant to be print
    l::Int -> Length of the strings to be print. If zero, prints max length
"""
function showdet(D::Determinant, l::Int = 0)

    if l == 0
        l = length(bitstring(D.α))
    end
    println("α: "*reverse(bitstring(D.α))[1:l])
    println("β: "*reverse(bitstring(D.β))[1:l])

end

end #Module
