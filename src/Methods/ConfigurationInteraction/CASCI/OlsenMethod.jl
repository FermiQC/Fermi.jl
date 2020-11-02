using TensorOperations
using Combinatorics
using SparseArrays
using ArnoldiMethod

function CASCI{T}(Alg::OlsenMethod) where T <: AbstractFloat

    @output "Getting molecule...\n"
    molecule = Molecule()
    @output "Computing AO Integrals...\n"
    #aoint = ConventionalAOIntegrals()

    @output "Calling RHF module...\n"
    refwfn = Fermi.HartreeFock.RHF(molecule)
    ints = refwfn.ints

    @output "Transforming Integrals for CAS computation...\n"
    # Read options
    frozen = Fermi.CurrentOptions["cas_frozen"]

    nmo = refwfn.ndocc + refwfn.nvir

    act_elec = 2*(refwfn.ndocc - frozen)

    if act_elec < 0
        error("\nInvalid number of frozen orbitals ($frozen) for $(2*refwfn.ndocc) electrons.")
    end

    # Active = -1 means FCI, with frozen
    if Fermi.CurrentOptions["cas_active"] == -1
        active = nmo - frozen
    else
        active = Fermi.CurrentOptions["cas_active"]
    end

    if active ≤ act_elec/2
        error("\nNumber of active orbitals ($active) too small for $(act_elec) active electrons")
    end

    if active+frozen > nmo
        error("\nNumber of active ($active) and frozen orbitals ($frozen) greater than number of orbitals ($nmo)")
    end

    s = 1:(frozen+active)
    h = T.(Fermi.Integrals.transform_fock(ints["T"] + ints["V"], ints.orbs["FU"][s], ints.orbs["FU"][s]))
    V = T.(Fermi.Integrals.transform_eri(ints["μ"], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s]))

    aoint = nothing
    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg)
end

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::OlsenMethod) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the Olsen Method.\n\n"


    nroot = Fermi.CurrentOptions["cas_nroot"]

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active
    
    @output "\nGenerating String replacement lists..."
    t = @elapsed begin
        strings, intree = get_αstrings(act_elec, active)
        Nstrings = length(strings)
    end
    
    @output " done in {} seconds.\n" t
    @output "\nNumber of Strings:      {:10d}\n" Nstrings
    @output "Number of Determinants: {:10d}\n" Nstrings^2

    @output "\nComputing auxiliar F matrix...\n"
    # Compute Fij
    F = h - 0.5*[tr(V[i,:,:,j]) for i=1:size(V,1), j=1:size(V,4)]
    @output "Computing Diagonal Hamiltonian matrix..."
    # Compute diagonal elements H(Iα, Iβ) 
    Hd = get_Hd(strings, h, V, act_elec)

    g = 2 # Number of initial guess vectors
    @output "\nUsing {:3d} unit vectors as initial guess" g
    b = zeros(Nstrings^2,g)

    for x=1:g
        b[x,x] = 1.0
    end

    #t = @elapsed begin
    #    Hfull = get_H_fromstrings(strings, intree, h, V, frozen)
    #    decomp, history = partialschur(Hfull, nev=nroot, tol=10^-12, which=LM())
    #    exact, ϕ = partialeigen(decomp)
    #end
    #@output "\n Exact FCI Energy: {:15.10f}\n" exact[1]+refwfn.molecule.Vnuc
    #@output "Time: {:5.5f}\n" t

    #b[1:100,1] .= ϕ[1:100,1]

    #t = @elapsed begin
    σ = zeros(size(b))
    σ = get_σ!(strings, intree, σ, b, F, V)
    #@assert 1==2
    #σ = Hfull*b
    G = b'*σ  # L x L
    λ = eigvals(G)[1] # Number
    α = eigvecs(G)[:,1] # L dimentional vector

    r = zeros(Nstrings^2,1)
    for i = 1:g
        r[:,1] += α[i]*(σ[:,i] - λ.*b[:,i])
    end
    δ = r ./(λ - sum(Hd)/length(Hd)) # N dimentional

    # Orthonormalize 
    b = gram_schmidt(b, δ)
    #gram_schmidt2(b, δ)

    ite = 1
    @output "\n\nIte   {:>15s}   {:>15s}   {:>15s}\n"  "Energy" "ΔE" "RMS"
    while ite < 50

        oldE = λ
        σ = zeros(size(b))
        σ = get_σ!(strings, intree, σ, b, F, V)
        #σ = Hfull*b
        G = b'*σ  # L x L
        λ = eigvals(G)[1] # Number

        ΔE = abs(oldE - λ)

        α = eigvecs(G)[:,1] # L dimentional vector

        r = zeros(Nstrings^2,1)
        for i = 1:g
            r[:,1] += α[i]*(σ[:,i] - λ.*b[:,i])
        end

        res = √(sum(r.^2))/Nstrings^2

        @output "{:3.0d}   {:15.10f}   {:15.10f}   {:15.10f}\n" ite λ ΔE res
        if res < 10^-10
            break
        end
        δ = r ./(λ .- Hd) # N dimentional
        #δ = r ./(λ - sum(Hd)/length(Hd)) # N dimentional

        # Orthonormalize 
        b = gram_schmidt(b, δ)
        #gram_schmidt2(b, δ)
        g = size(b,2)
        #if size(b,2) > 10
        #    old = deepcopy(b)
        #    b = zeros(Nstrings^2,1)
        #    for i = 1:g
        #        b[:,1] += α[i]*old[:,i]
        #    end
        #    g = 1
        #else
        #    g = size(b,2)
        #end
                

        ite += 1
    #end
    end
    @output "\n Davidson FCI Energy: {:15.10f}\n" λ+refwfn.molecule.Vnuc
    #@output "Time: {:5.5f}" t
end

function get_αstrings(Ne::Int, No::Int)                                             
                                                                                    
    Nae = Int(Ne/2)                                                                 
                                                                                    
    zeroth = vcat(repeat([true], Nae),repeat([false], No-Nae))                      
    conv = [2^i for i = 0:No-1]                                                     
                                                                                    
    perms = multiset_permutations(zeroth, length(zeroth))                           
                                                                                    
    string_list = [sum(conv[p]) for p in perms]                                     
                                                                                    
    # Initilize the interaction tree                                                
    intree = Array{Tuple{Int,Int,Int,Int},1}[]                                      
                                                                                    
    idx = 1                                                                         
    for p in perms                                                                  
        branch = Tuple{Int,Int,Int,Int}[]                                           
        for i = findall(p)                                                          
            for a = findall(.!p)                                                    
                                                                                    
                x,y = minmax(i,a)                                                   
                phase = (-1)^(sum(p[x:y])-1)                                        
                                                                                    
                newp = deepcopy(p)                                                  
                newp[i],newp[a] = false, true                                       
                newp = sum(conv[newp])                                              
                                                                                    
                address, = findall(x->x==newp, string_list)            
                push!(branch, (i,a,phase,address))                                  
            end                                                                     
            push!(branch, (i,i,1,idx))
        end                                                                         
        idx += 1
        push!(intree, branch)                                                       
    end                                                                             
                                                                                    
    return string_list, intree                                                      
end                                                                                 

function get_Hd(strings::Array{Int64,1}, h::Array{Float64,2}, V::Array{Float64,4}, Ne::Int) where N

    L = length(strings)
    Hd = zeros(L,L)
    αind = Array{Int64,1}(undef,Int(Ne/2))
    βind = Array{Int64,1}(undef,Int(Ne/2))

    for iα in 1:L
        orbindex!(strings[iα], αind)
        for iβ in 1:L
            orbindex!(strings[iβ], βind)
            Hd[iα, iβ] = Hd0(αind, βind, h, V)
        end
    end

    return reshape(Hd,L^2)
end

function orbindex!(D::Int, Out::Array{Int64,1})

    i = 1
    e = 1

    # Loop until 'e' electrons are found. Be careful! If the length of 'Out' is greater than
    # the number of electrons you will get stuck!
    while e ≤ length(Out)
        if 1<<(i-1) & D ≠ 0
            Out[e] = i
            e += 1
        end
        i += 1
    end
end

function get_σ!(strings::Array{Int64,1}, intree::Array{Array{NTuple{4,Int64},1},1}, σ::Array{Float64,2}, B::Array{Float64,2}, F::Array{Float64,2}, V::Array{Float64,4})

    L = length(strings)
    G = size(B,2)
    b = reshape(B, (L,L,G))
    σ = reshape(σ, (L,L,G))
    σ .= 0.0

    # Compute σ1
    for iα in 1:L
        # Kα is a single excitation from Iα 
        for (i,j,p1,kα) in intree[iα]
            # Jα is a single excitation from Kα
            elem = p1*F[i,j]
            for g=1:G, iβ=1:L
                σ[iβ,iα,g] += elem*b[iβ,kα,g]
                σ[iα,iβ,g] += elem*b[kα,iβ,g]
            end
            for (k,l,p2,jα) in intree[kα]
                elem = 0.5*p1*p2*V[i,j,k,l]
                for g=1:G, iβ=1:L
                    σ[iβ,iα,g] += elem*b[iβ,jα,g]
                    σ[iα,iβ,g] += elem*b[jα,iβ,g]
                end
            end
        end
    end

    # Compute σ3
    for iα in 1:L
        for (k,l,p1,jα) in intree[iα]
            for iβ in 1:L
                for (i,j,p2,jβ) in intree[iβ]
                    for g=1:G
                        σ[iβ,iα,g] += p1*p2*V[i,j,k,l]*b[jβ,jα,g]
                    end
                end
            end
        end
    end

    return reshape(σ, (L^2,G))
end

function gram_schmidt(b,δ)
    N = size(b,1)
    L = size(b,2)
    X = size(δ,2)

    δ = δ ./ norm(δ)

    out = Array{Float64,1}[]

    for i = 1:L
        push!(out, b[:,i])
    end

    for i = 1:X
        new = δ[:,i]
        for j = 1:length(out)
            new = new .- (new⋅out[j]).*out[j]
        end
        mag = norm(new)
        if mag > 10^-3
            push!(out, new/mag)
        end
    end

    return hcat(out...)
end

function gram_schmidt2(b,δ)
    N = size(b,1)
    L = size(b,2)
    X = size(δ,2)

    δ = δ ./ norm(δ)

    for i = 1:X
        new = δ[:,i]
        for j = 1:L
            new = new .- (new⋅b[:,j]).*b[:,j]
        end
        mag = norm(new)
        if mag > 10^-3
            b[:,i] += new/mag
        end
    end
end
