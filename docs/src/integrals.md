# Integrals and Orbitals

## IntegralHelper

The `Fermi.Integrals.IntegralHelper` object is designed to provide easy access to molecular 
integrals from Julia. 

### Creation
To create an `IntegralHelper`, you can simply call,

    Fermi.Integrals.IntegralHelper()

and one will be created based on the information in `Fermi.CurrentOptions`, which you may recall are set with `@set`. 
If you want more control, you can pass keyword arguments for the `Molecule` and `OrbDict` to be used,

    Fermi.Integrals.IntegralHelper(mol=mymol,orbs=myorbs)

### Use
Once an `IntegralHelper` object is created , integrals are simply accessed as,

    helper["S"] #AO basis overlap
    helper["T"] #AO basis kinetic
    helper["V"] #AO basis nuclear attraction
    helper["\mu"] #AO basis ERI
    helper["B"] #AO basis density fitted ERI (B(P,m,n))

    helper["OOOO"] # MO basis <ij|kl> TEI's 
    helper["OOVV"] # MO basis <ij|ab> TEI's
    helper["(OVOV)"] # MO basis (ia|jb) TEI's **chemists notation**
    helper["FOV"] # MO basis Fock matrix, OV block
    helper["BOV"] # MO basis DF TEI's, (ia|P)

Direct access to the underlying orbitals is possible,

    helper.orbs["O"]   # occupied orbitals, array of Orbital objects
    helper.orbs["[O]"] # occupied orbitals, collected into a matrix

To obtain unit normalized (non-CCA) integrals for d and higher angular momentum functions, call
`Fermi.Integrals.normalize(h::IntegralHelper,true)`. Note that _this will delete all current integrals_.

```@docs
Fermi.Integrals.IntegralHelper
```

```@docs
Fermi.Integrals.aokinetic
```

```@docs
Fermi.Integrals.aooverlap
```

```@docs
Fermi.Integrals.aodipole
```

```@docs
Fermi.Integrals.aonuclear
```

```@docs
Fermi.Integrals.aoeri
```

```@docs
Fermi.Integrals.dfaoeri
```

```@docs
Fermi.Integrals.aux_ri!
```

```@docs
Fermi.Integrals.aux_jk!
```

```@docs
Fermi.Integrals.getindex
```
