[![Build Status](https://travis-ci.com/FermiQC/Fermi.jl.svg?branch=master)](https://travis-ci.com/FermiQC/Fermi.jl)
[![Coverage Status](https://coveralls.io/repos/github/mdav2/Fermi.jl/badge.svg?branch=master)](https://coveralls.io/github/mdav2/Fermi.jl?branch=master)

# Fermi

Fermi is a quantum chemistry program written in (nearly) pure Julia. This code is developed at
the Center for Computational Quantum Chemistry at the University of Georgia. 

Fermi focuses on post Hartree--Fock methods. Currently, only restricted references are supported.
This is intended as a research code with an ever growing collection of methods implemented in
the package itself. However, the Fermi API is designed to make high performance pilot implementations
of methods achievable. 

Currently, we have implementations of:
- RHF (DF)
- RMP2 (DF)
- RMP3 (DF)
- RCCSD (DF)
- RCCSD(T)
- CASCI
- ACI
- ecCCSD
- ecCCSD(T) 

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
A minimal example of a computation is provided here. For more info check the doccumentation.

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
