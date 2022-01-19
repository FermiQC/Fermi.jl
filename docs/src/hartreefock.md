# Hartree-Fock 

The Hartree-Fock method is one of simplest *ab initio* methods where the wave function is modeled as a single Slater determinant of **spin-orbitals**

```math
\Psi_\text{HF}(1,2,...,N) = \frac{1}{\sqrt{N!}} \left| 
\begin{array}{c c c c}
\phi_1(1) & \phi_2(1) & ... & \phi_N(1) \\
\phi_1(2) & \phi_2(2) & ... & \phi_N(2) \\
\vdots & \vdots & \ddots & \vdots \\
\phi_1(N) & \phi_2(N) & ... & \phi_N(N) \\
\end{array}\right|
```
The spin-orbitals are constructed under the basis set approximation
```math
\phi_i(\vec{r}) = C_{\mu i} \chi_\mu(\vec{r})
```
> Note that sum over repeated indices is assumed.

where ``\chi(\vec{r})`` are contracted Gaussian basis functions from pre-constructed basis set such as `STO-3G`, `cc-pVDZ`, `ANO`, etc. 

These orbitals are constrained to be orthonormal to each other. Moreover, we choose to solve for the set of orbitals that diagonalize the Fock matrix
```math
D^\alpha_{\mu \nu} = C^\alpha_{\mu i} C^\alpha_{\nu i} \\[2mm]
F^\alpha_{\mu \nu} = H_{\mu \nu} + (D^\alpha_{\lambda\sigma} + D^\beta_{\lambda\sigma})(\mu \nu | \lambda \sigma) + D^\alpha_{\lambda\sigma}(\mu \nu | \lambda \sigma)
```
such that
```math
C^\alpha_{\mu i} F^\alpha_{\mu\nu} C^\alpha_{\nu j} = \delta_{ij} \epsilon_i
```
These are the canonical Hartree-Fock orbitals. 
> Equations also need to be solved for ``F^\beta``, but in the case of a restricted calculation, i.e. orbitals for both spins are taken to be the same, solving for ``\beta`` will yield the same results as for ``\alpha``.

## Restricted Hartree-Fock (RHF)

Minimal example
```julia
using Fermi

@molecule {
    He 0.0 0.0 0.0
}

@set basis 3-21g
@energy rhf
```
This will run a RHF computation on Helium using the `3-21G` basis set. Currently, Fermi does not support point group symmetry.

### Output file

The first part of the output gives an overview of the input information
```
He    0.000000000000    0.000000000000    0.000000000000


Charge: 0   Multiplicity: 1   
Nuclear repulsion:    0.0000000000
 Number of AOs:                            2
 Number of Doubly Occupied Orbitals:       1
 Number of Virtual Spatial Orbitals:       1
```

First the molecule XYZ is print. Followed by charge and multiplicity. Those will be taken as 0 and 1 by default, but can be controlled using `@set charge` and `@set multiplicity`. 

> ⚠️ RHF can only be used if the multiplicity is 1.

Next, we see the information about the iterations
```
 Iter.            E[RHF]         ΔE       Dᵣₘₛ        t     DIIS     damp
--------------------------------------------------------------------------------
    1     -2.8352184971  -2.835e+00   1.166e-01     0.78    false     4.71
    2     -2.8260289197   9.190e-03   2.885e-02     0.00    false     1.45
    3     -2.8157915919   1.024e-02   1.601e-02     0.00    false     0.00
    4     -2.8355956172  -1.980e-02   4.948e-02     0.18     true     0.00
    5     -2.8356798736  -8.426e-05   3.475e-03     0.00     true     0.00
    6     -2.8356798733   2.662e-10   8.346e-06     0.00     true     0.00
    7     -2.8356798736  -2.908e-10   6.418e-06     0.00     true     0.00
    8     -2.8356798736  -1.554e-14   4.527e-08     0.00     true     0.00
    9     -2.8356798728   8.546e-10   1.108e-05     0.14     true     0.00
    10    -2.8356798735  -7.070e-10   6.475e-06     0.00     true     0.00
    11    -2.8356798736  -1.477e-10   4.596e-06     0.00     true     0.00
    12    -2.8356798736  -8.882e-16   9.630e-09     0.00     true     0.00
    13    -2.8356798736   4.441e-16   1.687e-10     0.00     true     0.00
--------------------------------------------------------------------------------
 RHF done in  1.46s
```
Iterations are controlled using a few keywords. The convergence is achieved when
- The number of iterations reaches `scf_max_iter`

or
- ``\Delta E`` is less than `scf_e_conv` and ``D_{rms}`` is less than `scf_max_rms`

`DIIS` and `damp` are auxiliary strategies to reach convergency faster. 

Finally, the RHF energy is listed along with orbital energies
```
    @Final RHF Energy          -2.835679873641 Eₕ

   • Orbitals Summary

    Orbital            Energy    Occupancy
          1     -0.9035715084       ↿⇂
          2      2.0817026436         

   ✔  SCF Equations converged 😄
```

### RHF object

The computation returns a wave function object `Fermi.HartreeFock.RHF` which contains data useful for post-processing.

```@docs
Fermi.HartreeFock.RHF
```
