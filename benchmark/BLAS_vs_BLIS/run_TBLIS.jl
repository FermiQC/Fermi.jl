using Fermi

output("NUMBER OF THREADS: {}", Threads.nthreads())

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

MAX_CHAIN_LENGTH = 10
REPETITIONS = 10

timings = zeros(MAX_CHAIN_LENGTH, REPETITIONS)

@reset
@set {
    output TBLIS.out
    basis sto-3g
    tblis true
}

for N = 1:MAX_CHAIN_LENGTH

    get_hc(N)
    fc = 2*N
    Fermi.Options.set("drop_occ", fc)
    output("N = {:d}", N)

    # Cold run
    @set printstyle none
    ref = @energy rhf
    moints = Fermi.Integrals.IntegralHelper(orbitals=ref.orbitals)
    @energy moints=>CCSD

    # Replicatas
    @set printstyle file
    for l = 1:REPETITIONS
        timings[N,l] = @elapsed wf = @energy moints=>CCSD
        output("\n@@ TBLIS N = {:d}   l = {:d}  E = {:15.10f}   t = {:5.5f}\n", N, l, wf.energy, timings[N,l])
    end
end

output("\n Final Timings")
for N = 1:MAX_CHAIN_LENGTH
    t = sum(timings[N,:])/REPETITIONS
    output("N: {}   t(average) = {}", N, t)
end    
