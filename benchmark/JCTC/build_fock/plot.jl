using Plots
using Colors
using StatsPlots

include("collect_data.jl")

function plot_data()

    fcolor = ARGB(186/255,12/255,47/255, 0.8)
    pcolor = ARGB(39/255, 56/255, 150/255, 0.75)

    fermi = 1000*get_fermi()
    psi4 = 1000*get_psi4()

    nbas = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    bas = repr.(nbas)

    groupedbar([fermi psi4], c = [fcolor pcolor], legend=:topright, 
                label=["Fermi" "Psi4"], dpi=1200, orientation = :v, 
                bar_width=0.8,xticks=[2*i for i = 1:11], framestyle=:box)
    ylabel!("Time (milliseconds)")
    xlabel!("S22 Dimer")
end

function plot2_data()

    nbas = [58, 48, 104, 114, 264, 251, 321, 68, 96, 148, 228, 208, 264, 275, 321, 86, 138, 143, 147, 228, 275, 256]
    bas = repr.(nbas)

    alg_names = ["Fermi", "Psi4"]
    Nprogs = length(alg_names)
    prob_names = bas
    Nmols = length(prob_names)

    data = zeros(Nmols, Nprogs) 
    data[:,1] .= 1000*get_fermi()
    data[:,2] .= 1000*get_psi4()
    
    ctg = repeat(alg_names, inner = Nmols)
    nam = repeat(prob_names, outer = Nprogs)
    groupedbar(nam, data, group = ctg, xlabel = "S22 Dimer", ylabel = "Timing (milliseconds)",
        bar_width = 0.67,
        lw = 0, framestyle = :box, legend=:topright)
end
