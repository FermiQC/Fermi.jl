for N = 1:48
    run(`julia --threads $N run_X_threads.jl`)
end
