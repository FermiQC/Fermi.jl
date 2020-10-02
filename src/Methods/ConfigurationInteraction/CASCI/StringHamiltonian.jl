#using TensorOperations
#using Combinatorics
#using SparseArrays
#using ArnoldiMethod
#
#function CASCI{T}(Alg::SparseHamiltonian) where T <: AbstractFloat
#
#    @output "Getting molecule...\n"
#    molecule = Molecule()
#    @output "Computing AO Integrals...\n"
#    #aoint = ConventionalAOIntegrals()
#
#    @output "Calling RHF module...\n"
#    refwfn = Fermi.HartreeFock.RHF(molecule)
#    ints = refwfn.ints
#
#    @output "Transforming Integrals for CAS computation...\n"
#    # Read options
#    frozen = Fermi.CurrentOptions["cas_frozen"]
#
#    nmo = refwfn.ndocc + refwfn.nvir
#
#    act_elec = 2*(refwfn.ndocc - frozen)
#
#    if act_elec < 0
#        error("\nInvalid number of frozen orbitals ($frozen) for $(2*refwfn.ndocc) electrons.")
#    end
#
#    # Active = -1 means FCI, with frozen
#    if Fermi.CurrentOptions["cas_active"] == -1
#        active = nmo - frozen
#    else
#        active = Fermi.CurrentOptions["cas_active"]
#    end
#
#    if active ≤ act_elec/2
#        error("\nNumber of active orbitals ($active) too small for $(act_elec) active electrons")
#    end
#
#    if active+frozen > nmo
#        error("\nNumber of active ($active) and frozen orbitals ($frozen) greater than number of orbitals ($nmo)")
#    end
#
#    s = 1:(frozen+active)
#    h = T.(Fermi.Integrals.transform_fock(ints["T"] + ints["V"], ints.orbs["FU"][s], ints.orbs["FU"][s]))
#    V = T.(Fermi.Integrals.transform_eri(ints["μ"], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s], ints.orbs["FU"][s]))
#
#    aoint = nothing
#    CASCI{T}(refwfn, h, V, frozen, act_elec, active, Alg)
#end
#
#function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::SparseHamiltonian) where T <: AbstractFloat
#
#    # Print intro
#    Fermi.ConfigurationInteraction.print_header()
#    @output "\n    • Computing FCI with the SparseMatrix algorithm.\n\n"
#
#    nroot = Fermi.CurrentOptions["cas_nroot"]
#
#    @output "\n →  ACTIVE SPACE\n"
#    @output "Frozen Orbitals:  {:3d}\n" frozen
#    @output "Active Electrons: {:3d}\n" act_elec
#    @output "Active Orbitals:  {:3d}\n" active
#    
#    @output "\nGenerating α-strings.."
#    t = @elapsed begin
#        strings, tree = get_αstrings(act_elec, active)
#    end
#    
#    @output " done in {} seconds.\n" t
#    @output "\nNumber of α-strings:    {:10d}\n" length(strings)
#    @output "\nNumber of Determinants: {:10d}\n" length(strings)^2
#
#    @output "\nBuilding Sparse Hamiltonian..."
#end

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

function get_H_fromstrings(string_list::Array{Int64,1}, intree::Array{Array{NTuple{4,Int64},1},1}, h::Array{Float64,2}, V::Array{Float64,4}, frozen::Int)

    # Dense H for now
    L = length(string_list)
    H = zeros(L^2, L^2)

    # Build Aux matrix
    F = h - 0.5*[tr(V[i,:,:,j]) for i=1:size(V,1), j=1:size(V,4)]

    # Get 1-electron matrix elements
    for iα in 1:L
        for (i,j,p,jα) in intree[iα]
            elem = p*F[i,j]
            for iβ in 1:L
                d1 = L*(iα-1) + iβ
                d2 = L*(jα-1) + iβ
                H[d1,d2] += elem
            end
        end
    end

    for iβ in 1:L
        delem = 0.0
        for (i,j,p,jβ) in intree[iβ]

            elem = p*F[i,j]
            for iα in 1:L
                d1 = L*(iα-1) + iβ
                d2 = L*(iα-1) + jβ 
                H[d1,d2] += elem
            end
        end
    end

    #Get 2-electron matrix elements H[I,J] = 0.5*(ij,kl)*γijIK*γklKJ
    
    # αα excitations 
    for iα in 1:L
        d1 = L*(iα-1)
        # Kα is a single excitation from Iα 
        for (i,j,p1,kα) in intree[iα]
            # Jα is a single excitation from Kα
            for (k,l,p2,jα) in intree[kα]
                d2 = L*(jα-1)
                elem = 0.5*p1*p2*V[i,j,k,l]
                for iβ in 1:L
                    d1 += iβ
                    d2 += iβ
                    H[d1,d2] += elem
                    d1 -= iβ
                    d2 -= iβ
                end
            end
        end
    end

    # ββ excitations 
    for iβ in 1:L
        # Kβ is a single excitation from Iβ
        for (i,j,p1,kβ) in intree[iβ]
            # Jβ is a single excitation from Kβ
            for (k,l,p2,jβ) in intree[kβ]
                elem = 0.5*p1*p2*V[i,j,k,l]
                # Skip diagonal (I=J) cases for now
                #if jβ == iβ
                #    continue
                #end
                # I, K and J have the same alpha string
                for iα in 1:L
                    d = L*(iα-1)
                    d1 = d + iβ
                    d2 = d + jβ
                    H[d1,d2] += elem
                end
            end
        end
    end

    #αβ excitations
    for iα in 1:L
        # Kβ is equal Iβ    
        for iβ in 1:L
            # Kα is a single excitation from Iα and it is equal Jα
            for (i,j,p1,kα) in intree[iα]
                # Jβ is a single excitation from Kβ = Iβ
                for (k,l,p2,jβ) in intree[iβ]
                    elem = p1*p2*V[i,j,k,l]
                    d1 = L*(iα-1) + iβ
                    d2 = L*(kα-1) + jβ
                    H[d1,d2] += elem
                end
            end
        end

    end

    ##βα excitations
    #for iβ in 1:L
    #    # Kα is equal Iα
    #    for iα in 1:L
    #        # Kβ is a single excitation from Iβ
    #        for (i,j,p1,kβ) in intree[iβ]
    #            # Jα is a single excitation from Iα
    #            for (k,l,p2,jα) in intree[iα]
    #                elem = 0.5*p1*p2*V[i,j,k,l]
    #                d1 = L*(iα-1) + iβ
    #                d2 = L*(jα-1) + kβ
    #                H[d1,d2] += elem
    #            end
    #        end
    #    end
    #end

    return H
end
