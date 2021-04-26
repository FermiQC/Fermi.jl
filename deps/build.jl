# Fetch Libcint
run(`git clone https://github.com/sunqm/libcint.git`)
mkdir(joinpath(@__DIR__, "libcint/build"))
cd(joinpath(@__DIR__, "libcint/build"))

# Compile Libcint
println("Running CMAKE...")
run(`cmake ..`)
println("Running MAKE...")
run(`make`)

# Fetch binary 
println("Cleaning up...")
if Sys.islinux()
    run(`cp libcint.so ../../libcint.bin`)
elseif Sys.isapple()
    run(`cp libcint.dylib ../../libcint.bin`)
else
    error("Could not resolve OS")
end

# Untar basis set library and Clean up
cd(@__DIR__)
run(`tar -zxf lib.gz`)
rm(joinpath(@__DIR__, "libcint"), recursive=true, force=true)