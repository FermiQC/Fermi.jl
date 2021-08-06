<p align="center">
  <img src="docs/src/assets/logo.svg" width="400" alt=""/>
</p>

<table align="center">
  <tr>
    <th>Documentation</th>
    <th>Build Status</th>
    <th>License</th>
  </tr>
  <tr>
    <td align="center">
      <a href=https://FermiQC.github.io/Fermi.jl/dev>
      <img src=https://img.shields.io/badge/docs-dev-blue.svg>
      </a> 
    </td>
    <td align="center">
      <a href=https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml>
      <img src=https://github.com/FermiQC/Fermi.jl/actions/workflows/CI.yml/badge.svg>
      </a> 
      <a href=https://codecov.io/gh/FermiQC/Fermi.jl>
      <img src=https://codecov.io/gh/FermiQC/Fermi.jl/branch/master/graph/badge.svg?token=EWRG6Q7FK9>
      </a> 
    </td>
    <td align="center">
      <a href=https://github.com/FermiQC/Fermi.jl/blob/master/LICENSE>
      <img src=https://img.shields.io/badge/License-MIT-blue.svg>
      </a>
    </td>
  </tr>
</table>

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


## Contribute
PR's, issues, and suggestions are very welcome! You might consider reaching out before starting
work so that we can avoid duplication of efforts. Check the roadmap below for an idea of where this project is heading towards. Contact [Gustavo Aroeira](https://github.com/gustavojra) for any inquiries. 

## Roadmap

Fermi is a collection of *ab initio* methods. The long term goal is to provide production level implementations for daily applications.

### New features
1. Gradients for current methods. Three types of gradients are going to be considered: analytical, finite diferences, and automatic differentiation. 
2. Implementations for unrestricted refences, e.g. UCCSD(T).
3. High-order coupled cluster methods, such as CCSDT, and CCSDT(Q).
4. CBS extrapolation schemes and focal-point analysis.
5. Local correlation methods.
6. Interface with external modules for geometry optimization, vibrational analysis and themodynamic properties.

### Improvements
1. In an effort to improve the composability within the Julia chemistry community, some modules are going to be factorized out and Fermi will act as a high-level interface. For example, we intend to create a HartreeFock.jl package with all the current code in `Fermi.HartreeFock`. This way, anyone interest in the bare Hartree-Fock code can have a cleaner access to it. 

2. Performance boosts are always welcome! We need further testing and comparisons with well establish codes to find points to be improved. New backends for `BLAS` or `TBLIS` may also be considered.

## Employee of the month
   <img src="https://i.ibb.co/JQzmwTf/new.gif" alt="">
