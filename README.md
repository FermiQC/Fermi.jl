[![CI](https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FermiQC/Fermi.jl/branch/master/graph/badge.svg?token=EWRG6Q7FK9)](https://codecov.io/gh/FermiQC/Fermi.jl)

<img src="images/Fermi_logo.jpg" alt="Fermi Logo" width="350"/>

Fermi is a quantum chemistry program written in (nearly) pure Julia. This code is developed at
the [Center for Computational Quantum Chemistry](https://github.com/CCQC) at the University of Georgia under the supervision 
of [Dr. Justin M. Turney](https://github.com/jturney) and Prof. Henry F. Schaefer.

This work is supported by the U.S. National Science Foundation under grant number CHE-1661604.

### !! Warning
> Fermi is going through a major refactorization. Much of the code will be non functional. 
Please checkout at v0.2.0 version to try the code.

## Documentation

[DEV](https://FermiQC.github.io/Fermi.jl/dev)

Fermi focuses on post Hartree--Fock methods. Currently, only restricted references are supported.
This is intended as a research code with an ever growing collection of methods implemented in
the package itself. However, the Fermi API is designed to make high performance pilot implementations
of methods achievable. 

Currently, we have implementations of:

| Method    | Conv. | DF |
|-----------|-------|----|
| RHF       |  Y    |  Y |
<!--- | RMP2      |  Y    |  Y |
| RMP3      |  N    |  Y |
| RCCSD     |  Y    |  Y |
| RCCSD(T)  |  Y    |  Y |
| BCCD      |  Y    |  Y |
| CASCI     |  Y    |  N |
| ACI       |  Y    |  N |
| ecCCSD    |  Y    |  N |
| ecCCSD(T) |  Y    |  N |
-->



## Contribute
PR's, issues, and suggestions are very welcome! You might consider reaching out before starting
work so that we can avoid duplication of efforts.

## Install
Install Fermi by running,
```
pkg> add Fermi
```
If you would like the latest updates, use instead
```
pkg> add Fermi#master
```
Everything should work automatically on x86 architechtures for Linux and macOS. Windows is not
supported. You may run into issues when building Lints, the interface software between the
Libint2 integral code and Fermi. These errors can be a bit cryptic, so please reach out 
if you encounter any.

## Running single point energies
A minimal example of a computation is provided here. For more info check the documentation.

First, define a molecule
```
@molecule {
  O        1.2091536548      1.7664118189     -0.0171613972
  H        2.1984800075      1.7977100627      0.0121161719
  H        0.9197881882      2.4580185570      0.6297938830
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
