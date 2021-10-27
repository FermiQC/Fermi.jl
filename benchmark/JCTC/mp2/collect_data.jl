using PrettyTables

function raw_average(p = :fermi, engine = "MKL")
    if p == :fermi
        timings = get_fermi(engine)
    else
        timings = get_psi4()
    end
    avg = [sum(timings[s,:])/10 for s = 1:22]
    return avg
end

function get_fermi(engine = "MKL")
    timings = zeros(22,10)
    for s in 1:22
        for i in 1:10
            rex = Regex("$(engine) DF-MP2 Energy: S$s\\s+?$i run:\\s+([0-9]*\\.?[0-9]*)")
            timings[s,i] = get_linedata(rex, "fermi/output.dat")
        end
    end
    return timings
end

function get_psi4()
    Erex = r"DFMP2\sEnergy\s+:.+?([0-9]*\.?[0-9]*)w\s+[0-9]+\scalls"
    Rrex = r"DFMP2\sBia Read\s+:.+?([0-9]*\.?[0-9]*)w\s+[0-9]+\scalls"
    timings = zeros(22, 10)
    for s in 1:22
        for i in 1:10
            path = joinpath(@__DIR__, "psi4/R$i/S$s/timer.dat")
            timings[s,i] = get_linedata(Erex, path) - get_linedata(Rrex, path)
        end
    end
    return timings
end

function get_linedata(rex, fname)
    fpath = joinpath(@__DIR__, fname)
    out = 0.0
    f = false
    for l = eachline(fpath)
        m = match(rex, l)
        if m !== nothing
            out = m.captures[1] |> String |> x->parse(Float64,x)
            f = true
            break
        end
    end
    if !f
        error("Data not found for $(rex) in $(fname)")
    end
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    fermi_mkl = get_mkl()
    fermi_octavian = get_octavian()
    psi4 = get_psi4()
    
    @pt :header = ["Fermi - MKL", "Fermi - Octavian", "Psi4"] [fermi_mkl fermi_octavian psi4]
end
