using Fermi

function get_hc(N)

    molstring = ""
    f = false

    for l = eachline(joinpath(@__DIR__, "hydrocarbons.xyz"))
        if occursin("C$N", l)
            f = true
            continue   
        end

        if f && occursin("--", l)
            break
        end

        if f
            molstring *= l*"\n"
        end
    end

    Fermi.Options.set("molstring", molstring)
end

t_blas = []
e_blas = []
t_tblis = []
e_tblis = []

@reset
Nvals = collect(1:10)

@set {
    printstyle both
    output hydrocarbon_tblis.out
    basis cc-pvdz
}

for N = Nvals
    get_hc(N)
    output("N = {:d}", N)

    fc = 2*N

    Fermi.Options.set("drop_occ", fc)

    @set tblis false
    t = @elapsed e = @energy ccsd;
    push!(t_blas, t)
    push!(e_blas, e.energy)
    output("\n@@ BLAS N = {:d}   {:5.5f}   {:15.10f}\n", N, t, e.energy)

    @set tblis true
    t = @elapsed e = @energy ccsd;
    push!(t_tblis, t)
    push!(e_tblis, e.energy)
    output("\n@@ TBLIS N = {:d}   {:5.5f}   {:15.10f}\n", N, t, e.energy)
end