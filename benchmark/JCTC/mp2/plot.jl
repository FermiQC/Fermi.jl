using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots

include("collect_data.jl")


function plot_D()

    mcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)
    ocolor = :darkolivegreen
    

    mkl = get_mkl()
    oct = get_octavian()
    psi4  = get_psi4()

    nao = [144, 116, 236, 264, 592, 574, 724, 172, 232, 350, 528, 472, 592, 632, 724, 204, 322, 336, 338, 528, 632, 588]
    perm = sortperm(nao)
    nao = nao[perm]

    p1 = groupedbar([oct mkl psi4], c = [ocolor mcolor pcolor], legend=:topright, 
                label=["Fermi - Octavian" "Fermi - MKL" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,xticks=([2, 6, 12, 18, 22],["S2","S6", "S12","S18","S22"]), framestyle=:box)

    ylabel!("Time (seconds)")
    xlabel!("Dimer")
    p1
    p2 = scatter(nao, mkl[perm] ./ psi4[perm], label="", yticks=[0.5, 1.0, 1.5], ylims=[0.3,1.8], 
                color=:black, ylabel="Relative time\n(Fermi/Psi4)", framestyle=:box, marker=:diamond, markercolor=mcolor,
                xlabel="Number of basis functions", legend=:bottomleft)
    

    scatter!(p2, nao, oct[perm] ./ psi4[perm], label="", markercolor = ocolor, marker=:diamond)
    hline!(p2, [1], color=:gray, linestyle=:dot, label="")
    

    plot(p1, p2, layout = grid(2,1, heights=[0.7, 0.3]), leftmargin=0.8cm)
end

function plot_E()

    mcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)
    ocolor = :palegreen4
    

    mkl   = 1000 .* raw_average(:fermi, "MKL")
    oct   = 1000 .* raw_average(:fermi, "Octavian") 
    psi4  = 1000 .* raw_average(:psi4) 

    nao = [144, 116, 236, 264, 592, 574, 724, 172, 232, 350, 528, 472, 592, 632, 724, 204, 322, 336, 338, 528, 632, 588]
    perm = sortperm(nao)
    nao = nao[perm]

    p1 = groupedbar([oct mkl psi4], c = [ocolor mcolor pcolor], legend=:topright, 
                label=["Fermi - Octavian" "Fermi - MKL" "Psi4"], dpi=1200, orientation = :v, 
                yaxis=:log, 
                ylim = (10^-0.5,10^6),
                bar_width=0.8,xticks=([2, 6, 12, 18, 22],["S2","S6", "S12","S18","S22"]), framestyle=:box)

    ylabel!("Time (miliseconds)")
    xlabel!("Dimer")
    p1
    p2 = scatter(nao, mkl[perm] ./ psi4[perm], label="MKL", yticks=[0.5, 1.0, 1.5], ylims=[0.3,1.8], 
                color=:black, ylabel="Relative time\n(Fermi/Psi4)", framestyle=:box, marker=:diamond, markercolor=mcolor,
                xlabel="Number of basis functions", legend=:bottomleft, legendfontsize=6, legendfonthalign=:left)
    

    scatter!(p2, nao, oct[perm] ./ psi4[perm], label="Octavian", markercolor = ocolor, marker=:diamond)
    hline!(p2, [1], color=:gray, linestyle=:dot, label="")
    

    plot(p1, p2, layout = grid(2,1, heights=[0.7, 0.3]), leftmargin=0.8cm, dpi=1200)

    savefig("dfmp2.png")
end