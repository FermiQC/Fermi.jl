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

t_lints = []
e_lints = []
t_libcint = []
e_libcint = []

@reset
Nvals = collect(1:20)

# Cold run
@set printstyle none
@energy rmp2;
@set lints false
@energy rmp2;

@set {
    printstyle both
    df true
    basis cc-pvdz
    jkfit cc-pvtz-jkfit
    rifit cc-pvtz-rifit
}

for N = Nvals
    println(N)
    get_hc(N)

    @set lints true
    t = @elapsed e = @energy rmp2;
    push!(t_lints, t)
    push!(e_lints, e.energy)

    @set lints false
    t = @elapsed e = @energy rmp2;
    push!(t_libcint, t)
    push!(e_libcint, e.energy)
end

output("\n\n{:2s}   {:15s}   {:5s}   {:15s}   {:5s}   {:15s}   {:5s}", "N", "E(Lints)", "t(Lints)", "E(Libcint)", "t(libcint)", "ΔE", "Δt")
for i = eachindex(Nvals)

    output("{:2d}   {:15.10f}   {:5.5f}   {:15.10f}   {:5.5f}   {:15.10f}   {:5.5f}", 
    Nvals[i], e_lints[i], t_lints[i], e_libcint[i], t_libcint[i], e_libcint[i]-e_lints[i], t_libcint[i]-t_lints[i])
end
