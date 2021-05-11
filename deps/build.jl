# Fetch Libcint
run(`git clone https://github.com/sunqm/libcint.git`)
mkdir(joinpath(@__DIR__, "libcint/build"))
cd(joinpath(@__DIR__, "libcint/build"))

# Compile Libcint
@info " => Building Libcint"

@info " => Running CMAKE..."
run(`cmake ..`)
@info " => Running MAKE..."
run(`make`)

# Fetch binary 
@info " => Cleaning up..."
if Sys.islinux()
    run(`cp libcint.so ../../`)
elseif Sys.isapple()
    run(`cp libcint.dylib ../../`)
else
    @error "Could not resolve OS"
end

# Untar basis set library and Clean up
cd(@__DIR__)
run(`tar -zxf lib.gz`)
rm(joinpath(@__DIR__, "libcint"), recursive=true, force=true)