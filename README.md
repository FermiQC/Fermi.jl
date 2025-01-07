<p align="center">
  <img src="docs/src/assets/logo.svg" width="400" alt=""/>
</p>

<table align="center">
  <tr>
    <th>Documentation</th>
    <th>Build Status</th>
    <th>License</th>
    <th>Citation</th>
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
    <td align="center">
      <a href=https://pubs.acs.org/doi/10.1021/acs.jctc.1c00719>
      <img src=https://img.shields.io/badge/JCTC-10.1021/acs.jctc.1c00719-darkgreen.svg>
      </a>
    </td>
  </tr>
</table>

Fermi.jl is a quantum chemistry framework written in pure Julia. This code is developed at
the [Center for Computational Quantum Chemistry](https://github.com/CCQC) at the University of Georgia under the supervision 
of [Dr. Justin M. Turney](https://github.com/jturney) and Prof. Henry F. Schaefer.

This work is supported by the U.S. National Science Foundation under grant number CHE-1661604.

For an academic overview of this project, check our paper published in the Journal of Chemical Theory and Computation:
[Fermi.jl: A Modern Design for Quantum Chemistry](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00719)

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
| RCCSD(T)  |  Y    |  Y |
| FCI       |  Y    |  N |

## Install
Install Fermi by running,
```
pkg> add Fermi
```

> To access the package manager (`pkg>`) start the Julia terminal and hit `]`. 
> Alternatively, you can run
> ```
> julia> using Pkg
> julia> Pkg.add("Fermi")
> ```

If you would like the latest updates, use instead
```
pkg> add Fermi#master
```

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
Finally, run a computation
```
@energy ccsd;
```

## Fermi Ecosystem 

[Molecules.jl](https://github.com/FermiQC/Molecules.jl): A package that deals with `Atom` objects. It can parse XYZ files and compute properties related to the position of nuclei.

[GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl): A library for integrals and basis set objects. It can parse `.gbs` basis set files and create `BasisFunction` and `BasisSet` structures. Integrals over Gaussian basis are computed using `libcint`.

## Contribute
PR's, issues, and suggestions are very welcome! You might consider reaching out before starting
work so that we can avoid duplication of efforts. Besides browsing our issues for things that needs fixing/enhancement you can check the desired features below for ideas of what you can do! Contact [Gustavo Aroeira](https://github.com/gustavojra) for any inquiries. 

## Desired Features 

**1. New Methods for Energy Computations**
  * Unrestricted MP2
  * Unrestricted CC and High order CC methods (e.g. CCSDT, CCSDTQ)
  * Configuration Interaction (CID, CISD, FCI, selective CI, etc) 
  * Relativistic Methods
  * Methods beyond Born-Oppenheimer approximation
  * Excited States

**2. Properties**
  * Dipole Moments
  * Polarizability
  * IR intensities
  * Harmonic Frequencies
  * Anharmonic Frequencies

**3. Computational schemes**
  * GPU Computations
  * Distributed Computations
