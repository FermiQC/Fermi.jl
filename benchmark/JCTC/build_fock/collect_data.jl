using PrettyTables

function get_fermi()
    rex = r"@@@ Average Fock Time:\s([0-9]*\.?[0-9]*)"
    fpath = joinpath(@__DIR__, "fermi/output.dat")
    timings = zeros(22)
    i = 1
    for l = eachline(fpath)
        m = match(rex, l)
        if m !== nothing
            timings[i] = m.captures[1] |> String |> x->parse(Float64,x)
            i += 1
        end
    end

    @assert i == 23 # meaning that 22 entries were found

    return timings
end

function get_psi4()
    Frex = r"HF:\sForm\sF\s+:.+?([0-9]*\.?[0-9]*)w\s+([0-9]+)\scalls"
    Grex = r"HF:\sForm\sG\s+:.+?([0-9]*\.?[0-9]*)w\s+([0-9]+)\scalls"
    cd(joinpath(@__DIR__, "psi4")) 

    timings = zeros(22)

    # For each molecule
    for i = eachindex(timings)
        Ncalls = []
        Ft = []
        Gt = []
        for l = eachline("S$i/timer.dat")
            Fm = match(Frex, l)
            Gm = match(Grex, l)

            if Fm !== nothing
                t = Fm.captures[1] |> String |> x->parse(Float64,x)
                N = Fm.captures[2] |> String |> x->parse(Int,x)
                push!(Ncalls, N)
                push!(Ft, t)
            end

            if Gm !== nothing
                t = Gm.captures[1] |> String |> x->parse(Float64,x)
                N = Gm.captures[2] |> String |> x->parse(Int,x)
                push!(Ncalls, N)
                push!(Gt, t)
            end
        end
        # Make sure all data is consistent
        @assert all(Ncalls .== Ncalls[1])
        @assert all(Ft .== Ft[1])
        @assert all(Gt .== Gt[1])

        # Compute timing
        timings[i] = (Ft[1] + Gt[1]) / Ncalls[1]
    end

    cd("..")
    return timings
end

if abspath(PROGRAM_FILE) == @__FILE__
    fermi = get_fermi()
    psi4 = get_psi4()
    
    @pt :header = ["Fermi", "Psi4"] [fermi psi4]
end
