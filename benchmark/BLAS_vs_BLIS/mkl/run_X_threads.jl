using Fermi
using LinearAlgebra

Nt = Threads.nthreads()
Fermi.tblis_set_num_threads(Nt)
BLAS.set_num_threads(Nt)
output("NUMBER OF THREADS: {}", Threads.nthreads())

function get_hc(N)

    molstring = ""
    f = false

    for l = eachline(joinpath(@__DIR__, "../alkanes.xyz"))
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

MAX_CHAIN_LENGTH = 15

# Cold run
@set printstyle none
@energy CCSD

@reset
@set {
    basis sto-3g
    printstyle file
}

for N = 1:MAX_CHAIN_LENGTH

    get_hc(N)
    Fermi.Options.set("drop_occ", N)

    run(`echo TBLIS N = $N Threads = $Nt`)
    @set tblis true
    Fermi.Options.set("output", "TBLIS_C"*"$N"*"_T"*"$Nt"*".out")
    output("NUMBER OF THREADS: {}", Threads.nthreads())
    @energy ccsd

    run(`echo OpenBLAS N = $N Threads = $Nt`)
    @set tblis false
    output("NUMBER OF THREADS: {}", Threads.nthreads())
    Fermi.Options.set("output", "OPENBLAS_C"*"$N"*"_T"*"$Nt"*".out")
    @energy ccsd

end

using MKL

for N = 1:MAX_CHAIN_LENGTH

    get_hc(N)
    Fermi.Options.set("drop_occ", N)
    @set tblis false

    run(`echo MKL N = $N Threads = $Nt`)
    output("NUMBER OF THREADS: {}", Threads.nthreads())
    Fermi.Options.set("output", "MKL_C"*"$N"*"_T"*"$Nt"*".out")
    @energy ccsd
end
