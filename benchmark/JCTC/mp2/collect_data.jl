using PrettyTables

function get_fermi()
    rex = r"@@ Average time for DF-MP2 energy:\s*([0-9]*\.?[0-9]*)"
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
    Erex = r"DFMP2\sEnergy\s+:.+?([0-9]*\.?[0-9]*)w\s+[0-9]+\scalls"
    Rrex = r"DFMP2\sBia Read\s+:.+?([0-9]*\.?[0-9]*)w\s+[0-9]+\scalls"
    cd(joinpath(@__DIR__, "psi4")) 


    timings = zeros(22)
    for i = 1:10
        for s = 1:22
            for l = eachline("R$i/S$s/timer.dat")
                Em = match(Erex, l)
                Rm = match(Rrex, l)

                if Em !== nothing
                    timings[s] += Em.captures[1] |> String |> x->parse(Float64,x)
                elseif Rm !== nothing
                    timings[s] -= Rm.captures[1] |> String |> x->parse(Float64,x)
                    break
                end
            end
        end
    end
    return timings ./ 10
end

if abspath(PROGRAM_FILE) == @__FILE__
    fermi = get_fermi()
    psi4 = get_psi4()
    
    @pt :header = ["Fermi", "Psi4"] [fermi psi4]
end
