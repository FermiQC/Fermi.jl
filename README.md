# Fermi : Documentation
## Overview and Motivations

Fermi is a programming environment for writing arbitrary electronic structure and quantum chemical computations in the Julia programming language. Julia shows a lot of promise as a language for scientific computing, with many fields writing domain-specific applications in Julia. This project is intended to demonstrate some ways of working in this language, and showcase a proposed style of programming for expansion into a complete set of electronic structure programs.

In the benchmark and test folders, there are some files that act as examples for using the program. 

## Installing Fermi
These instructions are lifted from [this helpful site](https://tlienart.github.io/pub/julia/dev-pkg.html). 
In addition to the instructions below, there are some dependencies required. Please raise an issue if there is any difficulty with dependencies.
Make a directory where you will be placing Fermi. I'll use `<DEVDIR>` to represent that directory. Clone Fermi.


### TBLIS and Lints
Currently, Fermi depends on two unregistered packages TBLIS and Lints. Those can be added using their github address
```
julia> ]
(@v1.4) pkg> add https://github.com/FermiQC/Lints
(@v1.4) pkg> add https://github.com/FermiQC/TBLIS.jl
```
Make sure these packages are working by running
```
julia> ]
(@v1.4) pkg> test Lints
(@v1.4) pkg> test TBLIS
```
If an error like
```
ERROR: expected package `Lints [...]` to exist at path `...`
```
is raised, moving the package into `dev` may solve the problem
```
julia> ]
(@v1.4) pkg> dev Lints
```
### Making Julia aware of Fermi
First, clone Fermi
```
mkdir <DEVDIR>
cd <DEVDIR>
git clone https://github.com/mdav2/Fermi.jl.git
```
Now make the package manager Pkg aware of Fermi's presence.
```
julia> ]
(v1.3) pkg> dev <DEVDIR>/Fermi.jl
```
Now Fermi should be visible to Julia! To test this, `cd` and run,
```
$ julia
julia> using Fermi
```
You should see,
```
[ Info: Precompiling Fermi [9237668d-08c8-4784-b8dd-383aa52fcf74]
```
