# Backend

## Wavefunction.jl
`Wavefunction.jl` defines the `Wfn` type and constructors. These
constructors take in a Psi4 wavefunction object (PyObject) and convert
necessary information to Julia data structures.
```@docs
Fermi.Wavefunction
```

```@docs
Fermi.Wavefunction.Wfn(wfn)
```

```@docs
Fermi.Wavefunction.Wfn
```
## Transformation.jl

```@docs
Fermi.Transformation
```
## IntegralTransformation.jl
This module contain functions used to convert AO integrals into MO integrals.
It can be used to get ERI arrays and Fock matrices.
```@docs
Fermi.IntegralTransformation
```

```@docs
Fermi.IntegralTransformation.get_eri
```

```@docs
Fermi.IntegralTransformation.get_fock
```

```@docs
Fermi.Transformation.tei_transform
```

## DF.jl

```@docs
Fermi.DF
```
