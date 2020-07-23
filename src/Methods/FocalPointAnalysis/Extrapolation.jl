function extrapolate(data::Array{T,2}, basis::Array{String,1}, methods::Array{String,1}) where T <: AbstractFloat 
    X = []
    for b in basis
        if occursin("dz", lowercase(b))
            push!(X,2)
        elseif occursin("tz", lowercase(b))
            push!(X,3)
        elseif occursin("qz", lowercase(b))
            push!(X,4)
        else 
            m = match(r"(\d)z", lowercase(b))
            if m != nothing
                push!(X,parse(Int64,m.captures[1]))
            else
                error("Cardinality of $b could not be determined")
            end
        end
    end
    push!(X, Inf)

    # Create new data array
    nb = length(basis)
    nm = length(methods)
    extrapolated = zeros(nb+1, nm)
    extrapolated[1:nb, :] .= data
    extrapolated_indexes = Tuple{Int64,Int64}[]

    # Extrapolate SCF energy
    extrapolated[nb+1, 1], _, _ = feller_formula(data[nb-2:end,1], X[nb-2:end])
    push!(extrapolated_indexes, (nb+1,1))

    # Loop through correlated energies

    for m = 2:length(methods)
       
        corr = filter(x->x!=0, data[:, m])
        if length(corr) > 2
            # Subtract HF energy
            corr = corr .- data[1:length(corr), 1]    
            A,B = helgaker_formula(corr, X[1:length(corr)])

            for b = eachindex(X)
                if extrapolated[b,m] == 0.0
                    extrapolated[b,m] = A + B*X[b]^-3 + extrapolated[b,1]
                    push!(extrapolated_indexes,(b,m))
                end
            end
        else
            ΔE = extrapolated[1,m] - extrapolated[1,m-1]
            for b = 2:length(X)
                extrapolated[b,m] = extrapolated[b,m-1] + ΔE
                push!(extrapolated_indexes,(b,m))
            end
        end
    end

    return extrapolated, extrapolated_indexes
end

function feller_formula(E::Array{T,1}, X) where T <: AbstractFloat
    
    # SCF extrapolation: Escf = A + B*exp(-cX)
    # From: Feller 1993 (https://doi.org/10.1063/1.464749)

    # Get E1, E2, and E3
    l = length(E)
    E1 = E[l-2]
    E2 = E[l-1]
    E3 = E[l]

    X1 = X[l-2]
    X2 = X[l-1]
    X3 = X[l]

    C = log((E2 - E1)/(E3 - E2))
    B = (E3 - E2)/(exp(-C*X3) - exp(-C*X2))
    A = E3 - B*exp(-C*X3)
    
    return A,B,C
end

function helgaker_formula(E::Array{T,1}, X) where T <: AbstractFloat

    # Correlation energy extrapolation: Ecorr = A + B*X^-3
    # From Helgaker 1997 (https://doi.org/10.1063/1.473863)

    # Get E1, E2, and E3
    l = length(E)
    E1 = E[l-1]
    E2 = E[l]

    X1 = X[l-1]
    X2 = X[l]

    B = (E2 - E1)/(X2^-3 - X1^-3)
    A = E2 - B*X2^-3
    return A,B
end
