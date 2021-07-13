using Fermi

# Basic options trying to be uniform

Fermi.set_num_threads()

# Cold run
@set printstyle none
@energy mp2

function get_hc(N)

    molstring = ""
    f = false

    for l = eachline(joinpath(@__DIR__, "../alkanes/alkanes.xyz"))
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

@set {
    basis cc-pvdz
    scf_e_conv 1e-8
    scf_max_rms 1e-8
    cc_e_conv 1e-8
    cc_max_rms 1e-8
    tblis true
    printstyle file
    rifit cc-pvtz-rifit
    jkfit cc-pvtz-jkfit
}

prec = ["double", "single"]
dfs = [false, true]

for p in prec
    for d in dfs
        Fermi.Options.set("df", d)
        Fermi.Options.set("precision", p)
        Fermi.Options.set("output", p*"_"*repr(d)*".out")

        for i in 1:10

            # Read in xyz file
            get_hc(i)
            output("Alkane #{}", i)
        
            @freezecore
        
            # Run energy compt
            t = @elapsed begin
                @energy mp2
            end
            GC.gc()
            # Print out time
            output("@@@ Total Run Time: {}", t)
        end
    end
end
