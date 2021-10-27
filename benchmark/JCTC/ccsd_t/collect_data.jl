using PrettyTables
using Statistics

function hasoutliers(data)
    Qi, Qf = quantile(data, [0.01, 0.99])
    return any( data .> Qf) || any(data .< Qi)
end

function raw_average(p = :fermi, engine = "MKL")
    if p == :fermi
        timings = get_fermi(engine)
        avg = [sum(timings[s,:])/10 for s = 1:22]
    else
        timings = get_psi4()
        avg = [sum(timings[s,:])/5 for s = 1:22]
    end
    return avg
end

function get_fermi(engine = "MKL")
    timings = zeros(22,10)
    for s in 1:22
        for i in 1:10
            rex = Regex("S$s - \\(T\\) alg $(engine == "MKL" ? 1 : 2)\\s+?$i run:\\s+([0-9]*\\.?[0-9]*)")
            timings[s,i] = get_linedata(rex, "fermi/output.dat")
        end
    end
    return timings
end

function get_psi4()
    Erex = r"FNOCC:\striples\s+:.+?([0-9]*\.?[0-9]*)w\s+[0-9]+\scalls"
    timings = zeros(22, 5)
    for s in 1:22
        for i in 1:5
            path = joinpath(@__DIR__, "psi4/R$i/S$s/timer.dat")
            timings[s,i] = get_linedata(Erex, path)
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
