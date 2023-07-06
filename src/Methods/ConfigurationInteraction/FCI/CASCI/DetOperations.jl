module DetOperations

using Fermi
using LoopVectorization
using Combinatorics

export Determinant, get_determinants, αexcitation_level, βexcitation_level, phase, first_αexclusive, first_βexclusive
export αlist, βlist, αindex!, βindex!, second_αexclusive, second_βexclusive, detstring
export αocc!, βocc!, αvir!, βvir!, excitation_level, αexclusive, βexclusive, annihilate, create

"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant

structure holding alpha and beta strings for a determinant as integers.

# Fields

    α integer representing an alpha string
    β integer representing an beta string
"""
struct Determinant{T<:Integer}
    α::T
    β::T
end

"""
    Fermi.ConfigurationInteraction.DetOperations.Determinant(α::String, β::String)

Constructor function for Determinant object using strings.
"""
function Determinant(α::String, β::String; precision=Int64)
    
    αint = parse(precision, reverse(α); base=2) 
    βint = parse(precision, reverse(β); base=2) 
    return Determinant(αint, βint)
end

function get_determinants(Ne::Int, No::Int, nfrozen::Int)

    det_size = 
    if Fermi.Options.get("det_size") == 64
        Int64
    elseif Fermi.Options.get("det_size") == 128
        Int128
    else
        throw(ArgumentError("""Invalid determinant representation $(Fermi.Options.get("det_size"))"""))
    end

    Nae = Int(Ne/2)
    occ_string = repeat('1', nfrozen)

    zeroth = repeat('1', Nae)*repeat('0', No-Nae)

    perms = multiset_permutations(zeroth, length(zeroth))

    dets = [Determinant{det_size}[] for i = 1:Threads.nthreads()]
    @sync for αstring in perms
        Threads.@spawn begin
        for βstring in perms
            α = occ_string*join(αstring)
            β = occ_string*join(βstring)
            _det = Determinant(α, β;precision=det_size)
            push!(dets[Threads.threadid()], _det)
        end
        end #Threads.@spawn
    end

    dets = vcat(dets...)
    # Determinant list is sorted by its excitation level w.r.t the first determinant (normally HF)
    #sort!(dets, by=d->excitation_level(dets[1], d))

    return dets
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

    one = typeof(D.α)(1)
    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ length(Out)
        if one<<(i-1) & D.α ≠ 0
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

    one = typeof(D.α)(1)
    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If 'e' is greater than
    # the number of electrons you will get stuck!
    while e ≤ length(Out)
        if one<<(i-1) & D.β ≠ 0
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
    count_ones(D1.α ⊻ D2.α) >> 1
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexcitation_level(D1::Determinant, D2::Determinant)

Compare two determinants to return the excitation level of the beta electrons
"""
@inline function βexcitation_level(D1::Determinant, D2::Determinant;bail=false)
    count_ones(D1.β ⊻ D2.β) >> 1
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

    return αexclusive(D1.α, D2.α)
end

function αexclusive(D1α::Int64, D2α::Int64)

    αexcl = D1α ⊻ D2α & D1α

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

function αexclusive(D1α::Int128, D2α::Int128)

    αexcl = D1α ⊻ D2α & D1α

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while one(Int128)<<((i-1)&127) ≤ αexcl
        if one(Int128)<<((i-1)&127) & αexcl ≠ 0
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

    return first_αexclusive(D1.α, D2.α)
end

function first_αexclusive(D1α::Int64, D2α::Int64)

    x = D1α ⊻ D2α & D1α
    return 64 - leading_zeros(x&(~(x-1)))
end

function first_αexclusive(D1α::Integer, D2α::Integer)

    x = D1α ⊻ D2α & D1α
    return 128 - leading_zeros(x&(~(x-1)))
end

"""
    Fermi.ConfigurationInteraction.DetOperations.second_αexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the second alpha electron present in the second, but
not in the second
"""
function second_αexclusive(D1::Determinant, D2::Determinant)
    
    return second_αexclusive(D1.α, D2.α)
end

function second_αexclusive(D1α::Int64, D2α::Int64)

    x = D1α ⊻ D2α & D1α
    return 64 - leading_zeros(x)
end

function second_αexclusive(D1α::Integer, D2α::Integer)

    x = D1α ⊻ D2α & D1α
    return 128 - leading_zeros(x)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.βexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the indexes of the beta electrons present in the first, but
not in the second
"""
function βexclusive(D1::Determinant, D2::Determinant)

    return βexclusive(D1.β, D2.β)
end

function βexclusive(D1β::Int64, D2β::Int64)

    βexcl = D1β ⊻ D2β & D1β

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

function βexclusive(D1β::Int128, D2β::Int128)

    βexcl = D1β ⊻ D2β & D1β

    out = []
    i = 1
    # Save betas exclusives, in crescent order
    while one(Int128)<<(i-1) ≤ βexcl
        if one(Int128)<<(i-1) & βexcl ≠ 0
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

    return first_βexclusive(D1.β, D2.β)
end

function first_βexclusive(D1β::Int64, D2β::Int64)

    x = D1β ⊻ D2β & D1β
    return 64 - leading_zeros(x&(~(x-1)))
end

function first_βexclusive(D1β::Integer, D2β::Integer)

    x = D1β ⊻ D2β & D1β
    return 128 - leading_zeros(x&(~(x-1)))
end

"""
    Fermi.ConfigurationInteraction.DetOperations.second_βexclusive(D1::Determinant, D2::Determinant)

Compare two determinants to return the index of the second beta electron present in the second, but
not in the second
"""
function second_βexclusive(D1::Determinant, D2::Determinant)
    
    return second_βexclusive(D1.β, D2.β)
end

function second_βexclusive(D1β::Int64, D2β::Int64)

    x = D1β ⊻ D2β & D1β
    return 64 - leading_zeros(x)
end

function second_βexclusive(D1β::Integer, D2β::Integer)

    x = D1β ⊻ D2β & D1β
    return 128 - leading_zeros(x)
end

"""
    Fermi.ConfigurationInteraction.DetOperations.exclusive(D1::Determinant, D2::Determinant)

Returns a list with tuples (orbital index, spin) of orbitals populated in the first, but
not in the second Determinant.
"""
function exclusive(D1::Determinant, D2::Determinant)

    one = typeof(D1.α)(1)
    αexcl = D1.α ⊻ D2.α & D1.α
    βexcl = D1.β ⊻ D2.β & D1.β

    out = []
    i = 1
    # Save alphas exclusives, in crescent order
    while one<<((i-1)) ≤ αexcl
        if one<<((i-1)) & αexcl ≠ 0
            push!(out, (i, 'α'))
        end
        i += 1
    end
    i = 1
    # Save betas exclusives, in crescent order
    while one<<((i-1)) ≤ βexcl
        if one<<((i-1)) & βexcl ≠ 0
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

    one = typeof(D.α)(1)
    if spin == 'α'
        #if D.α & (1 << (orb-1)) == 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = 0
        i = typeof(D.α)(1)
        while i < (one << ((orb-1)))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α ⊻ (one << ((orb-1)))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        #if D.β & (1 << (orb-1)) == 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = typeof(D.β)(1)
        while i < (one << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β ⊻ (one << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

"""
    Fermi.ConfigurationInteraction.DetOperations.create(D::Determinant, orb::Int, spin::Char)

Creates a copy of the determinat where an electron was added to a specified orbital (orb) with a given spin (spin).
Returns the phase factor along with the new determinant
"""
function create(D::Determinant, orb::Int, spin::Char)

    one = typeof(D.α)(1)
    if spin == 'α'
        #if D.α & (1 << (orb-1)) ≠ 0
        #    error("Creation error. Orbital $orb is occupied")
        #end

        # Determine sign
        l = 0
        i = typeof(D.α)(1)
        while i < (one << (orb-1))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α | (one << (orb-1))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        #if D.β & (1 << (orb-1)) ≠ 0
        #    error("Annihilation error. Orbital $orb is not occupied")
        #end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = typeof(D.β)(1)
        while i < (one << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β | (one << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

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

    one = typeof(D.α)(1)
    e = 1
    for i in R
        if one<<(i-1) & D.α != 0
            Out[e] = i
            e += 1
        end
    end

end

function βocc!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    one = typeof(D.β)(1)
    e = 1
    for i in R
        if one<<(i-1) & D.β != 0
            Out[e] = i
            e += 1
        end
    end

end

function αvir!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    one = typeof(D.α)(1)
    e = 1
    for i in R
        if one<<(i-1) & D.α == 0
            Out[e] = i
            e += 1
        end
    end

end

function βvir!(D::Determinant, R::UnitRange{Int64}, Out::Array{Int64,1})

    e = 1
    one = typeof(D.β)(1)
    for i in R
        if one<<(i-1) & D.β == 0
            Out[e] = i
            e += 1
        end
    end
end

function detstring(D::Determinant,l::Int=64)

    αstr = ""
    βstr = ""
    str = ""

    i = 1
    for i = 1:l
        if 1<<(i-1) & (D.α & D.β) != 0
            αstr *= "1"
            βstr *= "1"
            str  *= "2"
        elseif 1<<(i-1) & D.α != 0
            αstr *= "1"
            βstr *= "0"
            str  *= "+"
        elseif 1<<(i-1) & D.β != 0
            αstr *= "0"
            βstr *= "1"
            str  *= "-"
        else
            αstr *= "0"
            βstr *= "0"
            str  *= "0"
        end
        i += 1
    end

    return "$str     $αstr     $βstr"
end
end #Module