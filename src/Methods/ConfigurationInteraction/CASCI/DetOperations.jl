using LoopVectorization
"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant(α::String, β::String)

Constructor function for Determinant object using strings.
"""
function Determinant(α::String, β::String)
    
    αint = parse(Int, reverse(α); base=2) 
    βint = parse(Int, reverse(β); base=2) 

    Determinant(αint, βint)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αlist(D::Determinant)

Return the alpha string of a Determinant as a list
"""
function αlist(D::Determinant)

    return [parse(Int, ss) for ss in reverse(bitstring(D.α))]
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βlist(D::Determinant)

Return the beta string of a Determinant as a list
"""
function βlist(D::Determinant)

    return [parse(Int, ss) for ss in reverse(bitstring(D.β))]
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αindexi!(D::Determinant, Out::Array{Int64,1})

Write the indexes of the occupid alpha electrons of the Determinant to a given list
"""
function αindex!(D::Determinant, Out::Array{Int64,1})

    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ length(Out)
        if 1<<(i-1) & D.α ≠ 0
            Out[e] = i
            e += 1
        end
        i += 1
    end
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βindex!(D::Determinant, Out::Array{Int64.1})

Write the indexes of the occupid beta electrons of the Determinant to a given list
"""
function βindex!(D::Determinant, Out::Array{Int64,1})

    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ length(Out)
        if 1<<(i-1) & D.β ≠ 0
            Out[e] = i
            e += 1
        end
        i += 1
    end
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αexcitation_level(D1::Determinant, D2::Determinant)

Compare two determinants to return the excitation level of the alpha electrons
"""
@inline function αexcitation_level(D1::Determinant, D2::Determinant;bail=false)
    count_ones(D1.α ⊻ D2.α)/2
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexcitation_level(D1::Determinant, D2::Determinant)

Compare two determinants to return the excitation level of the beta electrons
"""
@inline function βexcitation_level(D1::Determinant, D2::Determinant;bail=false)
    count_ones(D1.β ⊻ D2.β)/2
end

"""
    Fermi.ConfigurationInteraction.DetOperations.excitation_level(D1::Determinant, D2::Determinant)

Compare two determinants to return the excitation level between them
"""
@inline function excitation_level(D1::Determinant, D2::Determinant)
    return αexcitation_level(D1, D2) + βexcitation_level(D1,D2)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.αexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the indexes of the alpha electrons present in the first, but
not in the second
"""
function αexclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while 1<<((i-1)&63) ≤ αexcl
        if 1<<((i-1)&63) & αexcl ≠ 0
            push!(out, i)
        end
        i += 1
    end

    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.first_αexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the first alpha electron present in the first, but
not in the second
"""
function first_αexclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α

    i = 1
    # Save alphas exclusives, in crescent order
    while 1<<((i-1)) ≤ αexcl
        if 1<<((i-1)) & αexcl ≠ 0
            return i
        end
        i += 1
    end
    #return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.second_αexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the second alpha electron present in the second, but
not in the second
"""
function second_αexclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α

    i = 1
    sec = false
    # Save betas exclusives, in crescent order
    while 1<<(i-1) ≤ αexcl
        if 1<<(i-1) & αexcl ≠ 0
            if sec
                return i
            else
                sec = true
            end
        end
        i += 1
    end
    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the indexes of the beta electrons present in the first, but
not in the second
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
    Fermi.ConfigurationInteraction.DetOperations.first_βexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the first beta electron present in the first, but
not in the second
"""
function first_βexclusive(D1::Determinant, D2::Determinant)

    βexcl = D1.β ⊻ D2.β & D1.β

    i = 1
    # Save betas exclusives, in crescent order
    while 1<<((i-1)) ≤ βexcl
        if 1<<((i-1)) & βexcl ≠ 0
            return i
        end
        i += 1
    end

    #return out

end

"""
    Fermi.ConfigurationInteraction.DetOperations.second_βexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the second beta electron present in the second, but
not in the second
"""
function second_βexclusive(D1::Determinant, D2::Determinant)

    βexcl = D1.β ⊻ D2.β & D1.β

    i = 1
    sec = false
    # Save betas exclusives, in crescent order
    while 1<<((i-1)&63) ≤ βexcl
        if 1<<(i-1) & βexcl ≠ 0
            if sec
                return i
            else
                sec = true
            end
        end
        i += 1
    end
    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.exclusive(D1::Determinant, D2::Determinant)

Returns a list with tuples (orbital index, spin) of orbitals populated in the first, but
not in the second Determinant.
"""
function exclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α
    βexcl = D1.β ⊻ D2.β & D1.β

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while 1<<((i-1)) ≤ αexcl
        if 1<<((i-1)) & αexcl ≠ 0
            push!(out, (i, 'α'))
        end
        i += 1
    end
    i = 1
    # Save betas exclusives, in crescent order
    while 1<<((i-1)) ≤ βexcl
        if 1<<((i-1)) & βexcl ≠ 0
            push!(out, (i, 'β'))
        end
        i += 1
    end
    return out
end

"""
    Fermi.ConfigurationInteraction.DetOperations.annihilate(D::Determinant, orb::Int, spin::Char)

Creates a copy of the determinat where an electron from a specified orbital (orb) and spin (spin) was deleted. Returns the phase
factor along with the new determinant
"""
function annihilate(D::Determinant, orb::Int, spin::Char)

    if spin == 'α'
        #if D.α & (1 << (orb-1)) == 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = 0
        i = 1
        while i < (1 << ((orb-1)))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α ⊻ (1 << ((orb-1)))

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

Creates a copy of the determinat where an electron was added to a specified orbital (orb) with a given spin (spin).
Returns the phase factor along with the new determinant
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

Prints alpha and beta strings of a Determinant with a given length(l)
"""
function showdet(D::Determinant, l::Int = 0)

    if l == 0
        l = length(bitstring(D.α))
    end
    println("α: "*reverse(bitstring(D.α))[1:l])
    println("β: "*reverse(bitstring(D.β))[1:l])

end

function αocc!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    e = 1
    for i in R
        if 1<<(i-1) & D.α != 0
            Out[e] = i
            e += 1
        end
    end

end

function βocc!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    e = 1
    for i in R
        if 1<<(i-1) & D.β != 0
            Out[e] = i
            e += 1
        end
    end

end

function αvir!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    e = 1
    for i in R
        if 1<<(i-1) & D.α == 0
            Out[e] = i
            e += 1
        end
    end

end

function βvir!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    e = 1
    for i in R
        if 1<<(i-1) & D.β == 0
            Out[e] = i
            e += 1
        end
    end
end
