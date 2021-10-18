using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots

include("collect_data.jl")

function plot_A()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000*get_fermi()
    psi4  = 1000*get_psi4()

    groupedbar([fermi psi4], c = [fcolor pcolor], legend=:topright, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,xticks=[2*i for i = 1:11], framestyle=:box)
    ylabel!("Time (milliseconds)")
    xlabel!("S22 Dimer")

    savefig("A.png")
end

function plot_B()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000*get_fermi()
    psi4  = 1000*get_psi4()

    nao = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    perm = sortperm(nao)
    nao = nao[perm]
    fermi = fermi[perm]
    psi4 = psi4[perm]
    names = String[]

    for i = 1:22
        N = nao[i]
        S = perm[i]
        push!(names, "$S")
    end

    groupedbar(names, [fermi psi4], c = [fcolor pcolor], legend=:topleft, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,framestyle=:box,
                right_margin=1cm)

    ylabel!("Time (milliseconds)")
    xlabel!("S22 Dimer")
end

function plot_C()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000*get_fermi()
    psi4  = 1000*get_psi4()

    nao = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    perm = sortperm(nao)
    nao = nao[perm]
    fermi = fermi[perm]
    psi4 = psi4[perm]
    names = String[]
    snbf = String[]

    for i = 1:22
        N = nao[i]
        S = perm[i]
        push!(names, "S$S")
        push!(snbf, "$N")
    end

    p1 = groupedbar(names, [fermi psi4], c = [fcolor pcolor], legend=:topleft, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,framestyle=:box,
                right_margin=1cm)

    ylabel!("Time (milliseconds)")
    xlabel!("S22 Dimer")
    p2 = scatter(nao, fermi ./ psi4, label="", yticks=[0.5, 1.0, 1.5], ylims=[0.3,1.8], 
                color=:black, ylabel="Relative time\n(Fermi/Psi4)", framestyle=:box, marker=:utriangle, markercolor=:teal,
                xlabel="Number of basis functions")
    hline!(p2, [1], color=:gray, linestyle=:dot, label="")
    

    plot(p1, p2, layout = grid(2,1, heights=[0.7, 0.3]))
    #p1
end

function plot_D()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000*get_fermi()
    psi4  = 1000*get_psi4()

    nao = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    perm = sortperm(nao)
    nao = nao[perm]

    p1 = groupedbar([fermi psi4], c = [fcolor pcolor], legend=:topright, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,xticks=([2, 6, 12, 18, 22],["S2","S6", "S12","S18","S22"]), framestyle=:box)

    ylabel!("Time (milliseconds)")
    xlabel!("Dimer")

    p2 = scatter(nao, fermi[perm] ./ psi4[perm], label="", yticks=[0.5, 1.0, 1.5], ylims=[0.0,2.0], 
                color=:black, ylabel="Relative time\n(Fermi/Psi4)", framestyle=:box, marker=:diamond, markercolor=:navajowhite3,
                xlabel="Number of basis functions")
    hline!(p2, [1], color=:gray, linestyle=:dot, label="")
    

    plot(p1, p2, layout = grid(2,1, heights=[0.7, 0.3]), leftmargin=0.8cm)
end

function plot_E()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000 .* get_fermi()
    psi4  = 1000 .* get_psi4()

    nao = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    perm = sortperm(nao)
    nao = nao[perm]

    p1 = groupedbar(([fermi psi4]), c = [fcolor pcolor], legend=:topleft, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, yaxis=:log,
                ylim = (10^-0.5,10^4),
                bar_width=0.7,xticks=([2, 6, 12, 18, 22],["S2","S6", "S12","S18","S22"]), framestyle=:box)

    ylabel!("Time (milliseconds)")
    xlabel!("Dimer")

    p2 = scatter(nao, fermi[perm] ./ psi4[perm], label="", yticks=[0.5, 1.0, 1.5], ylims=[0.0,2.0], 
                color=:black, ylabel="Relative time\n(Fermi/Psi4)", framestyle=:box, marker=:diamond, markercolor=fcolor,
                xlabel="Number of basis functions")

    hline!(p2, [1], color=:gray, linestyle=:dot, label="")
    
    plot(p1, p2, layout = grid(2,1, heights=[0.7, 0.3]), leftmargin=0.8cm, dpi=1200)

    savefig("build_fock.png")
end