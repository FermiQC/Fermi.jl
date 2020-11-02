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

mutable struct SparseH
    ivals::Array{Int64,1}
    jvals::Array{Int64,1}
    vals::Array{Float64,1}
    size::Int64
end

function addtoindex!(SH::SparseH, val::Float64, i::Int, j::Int)

    if abs(val) < 10^-12
        return nothing
    end

    for x in 1:SH.size
        if (i,j) == (SH.ivals[x], SH.jvals[x])
            SH.vals[x] += val
            return nothing
        end
    end

    push!(SH.ivals, i)
    push!(SH.jvals, j)
    push!(SH.vals, val)
    SH.size += 1
end

function get_H_fromstrings(string_list::Array{Int64,1}, intree::Array{Array{NTuple{4,Int64},1},1}, h::Array{Float64,2}, V::Array{Float64,4}, frozen::Int)

    ## Dense H for now
    println("")
    L = length(string_list)

    @time begin
        H = zeros(L^2, L^2)
    end
    #H = SparseH([],[],[],0)
    #println("yay")

    ## Start with diagonal elements
    #@time begin
    #    ijvals = [(i,i) for i = 1:L^2]

    #    # Single excitations
    #    for s1 in 1:L
    #        for (i,j,p,excs1) in intree[s1]
    #            if i == j
    #                continue
    #            end
    #            for s2 in 1:L
    #                # α → α excitation
    #                d1 = L*(s1-1) + s2
    #                d2 = L*(excs1-1) + s2
    #                if !((d1,d2) in ijvals)
    #                    push!(ijvals, (d1,d2))
    #                end

    #                # β → β excitation
    #                d1 = L*(s2-1) + s1
    #                d2 = L*(s2-1) + excs1
    #                if !((d1,d2) in ijvals)
    #                    push!(ijvals, (d1,d2))
    #                end
    #            end
    #        end
    #    end

    #    # Double excitations (same spin)
    #    for s1 in 1:L
    #        for (i,j,p,excs1) in intree[s1]
    #            if i == j
    #                continue
    #            end
    #            for (k,l,p,dexcs1) in intree[s1]
    #                if k == l
    #                    continue
    #                end
    #                for s2 in 1: L
    #                    # αα → αα excitation
    #                    d1 = L*(s1-1) + s2
    #                    d2 = L*(dexcs1) + s2
    #                    if !((d1,d2) in ijvals)
    #                        push!(ijvals, (d1,d2))
    #                    end

    #                    # ββ → ββ excitation
    #                    d1 = L*(s2-1) + s1
    #                    d2 = L*(s2-1) + dexcs1
    #                    if !((d1,d2) in ijvals)
    #                        push!(ijvals, (d1,d2))
    #                    end
    #                end
    #            end
    #        end
    #    end

    #    # Double excitation (different spin)
    #    for iα in 1:L
    #        # Kβ is equal Iβ    
    #        for iβ in 1:L
    #            d1 = L*(iα-1) + iβ
    #            # Kα is a single excitation from Iα and it is equal Jα
    #            for (i,j,p1,kα) in intree[iα]
    #                if i == j
    #                    continue
    #                end
    #                # Jβ is a single excitation from Kβ = Iβ
    #                for (k,l,p2,jβ) in intree[iβ]
    #                if k == l
    #                    continue
    #                end
    #                    d2 = L*(kα-1) + jβ
    #                    if !((d1,d2) in ijvals)
    #                        push!(ijvals, (d1,d2))
    #                    end
    #                end
    #            end
    #        end
    #    end

    #    ivals = [ijvals[i][1] for i = 1:length(ijvals)]
    #    jvals = [ijvals[j][1] for j = 1:length(ijvals)]
    #    ijvals = 0.0
    #    vals = zeros(length(ivals))
    #    H = sparse(ivals, jvals, vals)
    #end

    # Build Aux matrix
    F = h - 0.5*[tr(V[i,:,:,j]) for i=1:size(V,1), j=1:size(V,4)]


    @time begin 
    # Get 1-electron matrix elements
    for iα in 1:L
        _d1α = L*(iα-1)
        for (i,j,p,jα) in intree[iα]
            _d2α = L*(jα-1)
            elem = p*F[i,j]
            for iβ in 1:L
                _xβ = L*(iβ-1)

                d1α = _d1α + iβ
                d2α = _d2α + iβ
                d1β = _xβ + iα
                d2β = _xβ + jα

                H[d1α,d2α] += elem
                H[d1β,d2β] += elem
                #addtoindex!(H, elem, d1α, d2α)
                #addtoindex!(H, elem, d1β, d2β)
                #push!(ivals, d1α)
                #push!(jvals, d2α)
                #push!(vals, elem)
                #push!(ivals, d1β)
                #push!(jvals, d2β)
                #push!(vals, elem)
            end
        end
    end
    end

    #Get 2-electron matrix elements H[I,J] = 0.5*(ij,kl)*γijIK*γklKJ
    @time begin 
    # αα excitations 
    for iα in 1:L
        _d1α = L*(iα-1)
        # Kα is a single excitation from Iα 
        for (i,j,p1,kα) in intree[iα]
            # Jα is a single excitation from Kα
            for (k,l,p2,jα) in intree[kα]
                _d2α = L*(jα-1)
                elem = 0.5*p1*p2*V[i,j,k,l]
                for iβ in 1:L

                    d1α = _d1α + iβ
                    d2α = _d2α + iβ

                    H[d1α,d2α] += elem
                    #addtoindex!(H, elem, d1α, d2α)

                    _xβ = L*(iβ-1)
                    d1β = _xβ + iα
                    d2β = _xβ + jα

                    H[d1β,d2β] += elem
                    #addtoindex!(H, elem, d1β, d2β)

                    #push!(ivals, d1α)
                    #push!(jvals, d2α)
                    #push!(vals, elem)

                    #push!(ivals, d1β)
                    #push!(jvals, d2β)
                    #push!(vals, elem)
                end
            end
        end
    end
    end

    #αβ excitations (ij|kl) <Iα|Eij|Kα> <Kβ|Ekl|Jβ> δ(Iβ|Kβ) δ(Kα|Jα)
    @time begin
    for iα in 1:L
        # Kβ is equal Iβ    
        for iβ in 1:L
            d1 = L*(iα-1) + iβ
            # Kα is a single excitation from Iα and it is equal Jα
            for (i,j,p1,kα) in intree[iα]
                # Jβ is a single excitation from Kβ = Iβ
                for (k,l,p2,jβ) in intree[iβ]
                    elem = p1*p2*V[i,j,k,l]
                    d2 = L*(kα-1) + jβ
                    H[d1,d2] += elem
                    #addtoindex!(H, elem, d1, d2)

                    #push!(ivals, d1)
                    #push!(jvals, d2)
                    #push!(vals, elem)
                end
            end
        end
    end
    end

    return sparse(H)
end
