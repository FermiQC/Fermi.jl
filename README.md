<p align="center">
  <img src="images/Fermi_logo.png" />
</p>

| **Documentation**| **Build** | **License**|
|:--:|:--:|:--:| 
|[![DEV][ddocs-img]][ddocs-url] | [![CI][ci-img]][ci-url] [![codecov][ccov-img]][ccov-url]| [![MIT license][lic-img]][lic-url]

[ddocs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[ddocs-url]: https://FermiQC.github.io/Fermi.jl/dev

[ci-img]: https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml

[ccov-img]: https://codecov.io/gh/FermiQC/Fermi.jl/branch/master/graph/badge.svg?token=EWRG6Q7FK9
[ccov-url]: https://codecov.io/gh/FermiQC/Fermi.jl

[lic-img]:https://img.shields.io/badge/License-MIT-blue.svg
[lic-url]:https://github.com/FermiQC/Fermi.jl/blob/master/LICENSE

Fermi is a quantum chemistry program written in (nearly) pure Julia. This code is developed at
the [Center for Computational Quantum Chemistry](https://github.com/CCQC) at the University of Georgia under the supervision 
of [Dr. Justin M. Turney](https://github.com/jturney) and Prof. Henry F. Schaefer.

This work is supported by the U.S. National Science Foundation under grant number CHE-1661604.

Fermi focuses on post Hartree--Fock methods. Currently, only restricted references are supported.
This is intended as a research code with an ever growing collection of methods implemented in
the package itself. However, the Fermi API is designed to make high performance pilot implementations
of methods achievable. 

Currently, we have implementations of:

| Method    | Conv. | DF |
|-----------|-------|----|
| RHF       |  Y    |  Y |
| RMP2      |  Y    |  Y |
| RCCSD     |  Y    |  Y |
| RCCSD(T)  |  Y    |  Y |


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
pkg> add Fermi
```
Everything should work automatically, the most flagile part is building the integral library [libcint](https://github.com/sunqm/libcint). The file [`deps/build.jl`](https://github.com/gustavojra/Fermi.jl/blob/master/deps/build.jl) contains simple commands to clone and build this library, you might need to modify it to better suit your system. If you do, rerun the build step using `pkg> build Fermi`. Please reach out 
if you encounter any problem.

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
