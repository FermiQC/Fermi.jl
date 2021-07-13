using Fermi
using LinearAlgebra

Nt = Threads.nthreads()
Fermi.tblis_set_num_threads(Nt)
BLAS.set_num_threads(Nt)

function get_hc(N)

    molstring = ""
    f = false

    for l = eachline(joinpath(@__DIR__, "alkanes.xyz"))
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

MAX_CHAIN_LENGTH = 20

# Cold run
@set printstyle none
@energy rhf

@reset
@set {
    basis cc-pvtz
    scf_e_conv 1e-8
    scf_max_rms 1e-8
    printstyle file
    output output.dat
}

for N = 1:MAX_CHAIN_LENGTH

    get_hc(N)
    # Run energy compt
    t = @elapsed begin
        @energy rhf
    end
    GC.gc()
    # Print out time
    output("@@@ Total Run Time: {}", t)
end
