# Fermi
## Overview and Motivations

Fermi is a programming environment for writing arbitrary electronic structure and quantum chemical computations in the Julia programming language. Julia shows a lot of promise as a language for scientific computing, with many fields writing domain-specific applications in Julia. This project is intended to demonstrate some ways of working in this language, and showcase a proposed style of programming for expansion into a complete set of electronic structure programs.

In the benchmark and test folders, there are some files that act as examples for using the program. 

## Installing Fermi
These instructions are lifted from [this helpful site](https://tlienart.github.io/pub/julia/dev-pkg.html). 
Currently, Fermi is unregistered, thus we use it with the `dev` mode.
In addition to the instructions below, there are some dependencies required. Please raise an issue if there is any difficulty with dependencies.
Make a directory where you will be placing Fermi. I'll use `<DEVDIR>` to represent that directory. Clone Fermi.

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

## Running single point energies
A minimal example of a computation is provided here. For more info check the doccumentation (TODO)
First, define a molecule
```
@molecule {
  O        1.2091536548      1.7664118189     -0.0171613972
  H        2.1984800075      1.7977100627      0.0121161719
  H        0.9197881882      2.4580185570      0.629793883
}
```
Choose a basis set
```
@set basis sto-3g
```
Finally run a computation
```
@energy ccsd;
```
