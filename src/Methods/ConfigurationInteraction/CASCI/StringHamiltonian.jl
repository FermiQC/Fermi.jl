using TensorOperations
using Combinatorics
using SparseArrays
using ArnoldiMethod

function CASCI{T}(Alg::SparseHamiltonian) where T <: AbstractFloat

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

function CASCI{T}(refwfn::Fermi.HartreeFock.RHF, h::Array{T,2}, V::Array{T,4}, frozen::Int, act_elec::Int, active::Int, Alg::SparseHamiltonian) where T <: AbstractFloat

    # Print intro
    Fermi.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the SparseMatrix algorithm.\n\n"

    nroot = Fermi.CurrentOptions["cas_nroot"]

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active
    
    @output "\nGenerating α-strings.."
    t = @elapsed begin
        dets = get_determinants(act_elec, active, frozen)
        Ndets = length(dets)
    end
    
    @output " done in {} seconds.\n" t
    @output "\nNumber of Determinants: {:10d}\n" Ndets

    @output "\nBuilding Sparse Hamiltonian..."
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
        end                                                                         
        push!(intree, branch)                                                       
    end                                                                             
                                                                                    
    return string_list, intree                                                      
end                                                                                 
