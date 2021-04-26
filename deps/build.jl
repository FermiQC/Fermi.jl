# Fetch Libcint
run(`git clone https://github.com/sunqm/libcint.git`)
mkdir(joinpath(@__DIR__, "libcint/build"))
cd(joinpath(@__DIR__, "libcint/build"))

# Compile Libcint
println("Running CMAKE...")
run(`cmake ..`)
println("Running MAKE...")
run(`make`)

# Clean up and untar basis set files
println("Cleaning up...")
run(`cp libcint.so ../../`)
cd(@__DIR__)
run(`tar -zxf lib.gz`)
rm(joinpath(@__DIR__, "libcint"), recursive=true, force=true)
rm(joinpath(@__DIR__, "lib.gz"), recursive=true, force=true)