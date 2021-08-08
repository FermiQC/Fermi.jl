for N = 2:48
    run(`julia --threads $N run_X_threads.jl`)
    run(`mkdir N$N`)

    for f in readdir(@__DIR__)
        if occursin(".out", f)
            mv(f, "N$N/$f")
        end
    end
end
