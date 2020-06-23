# Fermi : Documentation
## Overview and Motivations
> **(!)** This project is preliminary. I'm still quite new to Julia, and any feedback is very welcome! 

Fermi (pronounced "juice") is a programming environment for writing arbitrary electronic structure and quantum chemical computations in the Julia programming language. Julia shows a lot of promise as a language for scientific computing, with many fields writing domain-specific applications in Julia. This project is intended to demonstrate some ways of working in this language, and showcase a proposed style of programming for expansion into a complete set of electronic structure programs.

In the benchmark and test folders, there are some files that act as examples for using the program. 

## Installing Fermi
These instructions are lifted from [this helpful site](https://tlienart.github.io/pub/julia/dev-pkg.html). 
In addition to the instructions below, there are some dependencies required. Please raise an issue if there is any difficulty with dependencies.
Make a directory where you will be placing Fermi. I'll use `<DEVDIR>` to represent that directory. Clone Fermi.

### Linking to Psi4Numpy
> **(!) Temporary** Due to a conflict between Numpy and MKL, if you want to use Intel MKL, which can boost performance for tensor contractions, currently you must obtain MKL.jl from [this repo](https://github.com/fgerick/MKL.jl) NOT the JuliaComputing master branch (currently). 

Please do the following to make the Psi4Numpy interface visible to Julia.
```
conda create -n p4env python=3.7 psi4 -c psi4
```
To get the path to the Python executable,
```
$ conda activate p4env
$ which python
```
Then make this python visible to Julia:
```julia-repl
julia>ENV["PYTHON"] = <path-to-p4env-python>
julia>] build PyCall
```

### Making Julia aware of Fermi
```
mkdir <DEVDIR>
cd <DEVDIR>
git clone https://github.com/mdav2/Fermi.jl.git
```
Now make the package manager Pkg aware of Fermi' presence.
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

I strongly encourage running the test suite immediately! Pkg once again makes this easy,
```
julia> ]
(v1.3) pkg> test Fermi
```
or if you prefer
```
julia> import Pkg
julia> Pkg.test("Fermi")
```
You will see a bunch of package versions spill out, and the REPL will hang while it runs tests. This takes about 2 minutes on a very slow machine (Intel Core m3 processor). Finally you will see an output from the test suite. If all is well, you will see something like,
```
Test Summary:             | Pass  Fail  Error  Total
Fermi                      |   18                  18
```
The number of tests may be different, though. However, in the more likely case that something has an unfixed bug, in this case an issue with DiskTensors and an issue with MP2. "Fail" means the test ran without errors (syntax etc.) but did not produce the prescribed output. Usually this is a problem with the implementation of the method - please place an issue on GitHub! "Error" means the mistake was more boneheaded - usually because a breaking change was made without correctly modifying the test suite. 
```
Test Summary:             | Pass  Fail  Error  Total
Fermi                      |   16     1      1     18
  Wavefunction            |    4                   4
  CISingles               |    1                   1
  DiskTensors             |    6            1      7
    Smoke                 |    3                   3
    Dot                   |    3            1      4
      DiskVector          |    1            1      2
      DiskMatrix          |    1                   1
      DiskFourTensor      |    1                   1
  CoupledCluster          |    1                   1
  Integral Transformation |    2                   2
  MP2                     |    2     1             3
```

## Design
### Integral backends
This project uses other electronic structure programs to compute basic quantities like the one- and two-electron Hamiltonian integrals, as well as their counterparts for other operators. There is currently no intention of writing an integrals code specifically for Fermi.

As of now there is an interface to the Psi4 programs via the `psi4numpy` interface. This is simply the interface that I know the best, and should be extended in the future. There are plans for interfacing with the PySCF project and the NWChem project. If someone with knowledge wants to implement an interface to other programs e.g. Q-Chem, CFOUR, ORCA, Turbomole, or MOLPRO, they are very welcome. It is intended that the most robust interfaces should be to free and/or open source programs. 

Through the `psi4numpy` interface, RHF and UHF wavefunctions are obtained. For cases
of problematic spin contamination, the CUHF reference may be used to form a spin-restricted
UHF solution. 

### Naming
It is customary for Julia modules and data structures to have CamelCase names, such as `CoupledCluster.jl`. Please follow this aesthetically pleasing convention! 

## Submodules
This section contains a description of some preliminary modules in the Fermi environment. Many aspects are aspirational, and the description is not so much a description of current functionality as a statement of intent.
### Wavefunction.jl
Wavefunction.jl is the foundational module of Fermi. All Fermi programs will make use of Wavefunction.jl at some point. This module is the point of interaction between integral backends (e.g. `psi4numpy`) and the Fermi programming environment. An interface to an integral backend should produce a complete Wfn structure from the relevant sources. A Wfn structure is a representation of a reference determinant, such as a set of Hartree-Fock or DFT orbitals. 
> **Help** This module is functioning, and has a test set. Fleshing out documentation and contributing tests would be a great help, but other modules are in greater need.
### Determinant.jl
General representation of a Slater determinant. Used in modules CISingles.jl. 
> **(!) Help** This module is severely lacking - excellent starting project for an experienced programmer.
### DiskTensors.jl
This module describes a way of storing vectors, matrices, and rank four tensors on disk in a convenient way for use in electronic structure computations. 
This module uses the HDF5 binary format for storing and accessing arrays. This was chosen for its convenient interface, support for compressed I/O, and good support in Julia. 
> **Help** This module has all intended types (rank 1,2,4 tensors) implemented but many operations (+,-,/,',...) are not defined. Test suite is reasonably complete and documentation is rough.
### Davidson.jl
This module implements a simple Davidson solver, which currently has some unidentified bug. Use the IterativeSolvers.jl LOBPCG routine for in-core computations.
> **(!) Help** Rewriting the Davidson code is probably a good idea. A routine to collapse the trial vector subspace would make this module much more functional. A generalized implementation for non-symmetric matrices is required before EOM codes can be useful. 
### Direct.jl
This module contains necessary code for integral direct computations. Currently only interfaces to the Psi4 programs, but additions of interfaces are welcome. 
> **(!) Help** Function for contracting integrals for a fixed index is implemented. For most purpses, storing a 3-index quantity in memory is feasible; so incorporating this into energy routines would be very helpful.
### MatrixElement.jl
This module defines an interface for obtaining matrix elements for CI matrices.
>**(!) Help** This is just a skeleton at this point. Contributions to this module will greatly help a functioning FCI and arbitrary order CI code. Basic equations and citation to Szabo and Ostlund are in docstrings to help the intrepid contributor. 
### MollerPlesset.jl
Routines for Moller-Plesset perturbation theory computations are implemented here. Currently in-core and disk based RMP2 and UMP2 are implemented. 
>**(!) Help** Direct UMP2 would be an excellent contribution. 
### CISingles.jl
Specialized routines for computing configuration-interaction singles excited state wavefunctions are defined here. Corrections such as CIS(D) and variants defined here as well. Keep seperate from general CI code. Only in-core RCIS is implemented, and is currently not working.
>**(!) Help** UCIS, disk-based, and direct implementations are excellent targets. 
### CoupledCluster.jl
Routines for computing ground state coupled cluster energies are contained here. Currently there is RHF-CCD and RHF-CCSD implemented. UHF and ROHF codes are in the works.
>**(!) Help** Working on UHF-CC codes would be greatly appreciated. Coding a perturbative triples correction would be a straightforward addition. CCD and CCSD code could be updated to be able to use DiskTensors for Fxy and Wvxyz intermediates.
