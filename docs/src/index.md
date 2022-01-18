# Fermi.jl
Fermi is Julia framework for *ab initio* quantum chemistry. The main goal of Fermi is to offer two key resources:

- A collection of standard methods in electronic structure theory, such as Hartree--Fock and Coupled Cluster.
- An efficient development platform for quantum chemistry methods.

📝 For an academic overview of the project, check our paper published at the Journal of Chemical Theory and Computation:

[Fermi.jl: A Modern Design for Quantum Chemistry](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c00719)

# Installation

Fermi is a registered Julia package and, as such, it can be obtained using the standard Julia package manager. From the Julia terminal, use the `]` to move to the pkg manager
```julia
julia> # This is the standard Julia terminal, hit ] to go into Pkg
(@v1.6) pkg> # This is the package manager! Hit back space to leave this mode
```
Next add Fermi to the current environment.
```julia
(@v1.6) pkg> add Fermi
```
All the dependencies are going to be downloaded and installed and the code should be ready to work. To test the package you can run
```julia
(@v1.6) pkg> test Fermi
```

If you want a version of Fermi that you can modify, clone and check it for development
```julia
shell> git clone https://github.com/FermiQC/Fermi.jl
shell> cd Fermi.jl
(@v1.6) pkg> dev .
```

# Usage

Fermi can be used interactively through the Julia terminal, or you can write a Julia script which will act as the traditional input file present in other quantum chemistry packages. For example, a minimal script to run a RHF computation on a water molecule can be written as
```julia
using Fermi

@molecule {
  O        1.2091536548      1.7664118189     -0.0171613972
  H        2.1984800075      1.7977100627      0.0121161719
  H        0.9197881882      2.4580185570      0.6297938830
}

@set basis sto-3g
@energy rhf
```
If you save this file as `input.jl`, you can run it as a regular script
```
shell> julia --threads N input.jl
``` 
where N is the desired number of threads. Alternatively, you can set `export JULIA_NUM_THREADS=N` in your path.

Multiple options can be set using the macro `@set` with curly braces
```julia
@set {
  basis cc-pvdz
  df true
  jkfit cc-pvtz-jkfit
}
```

The macro `@energy` returns a wave function object associated with the methods requested. See the documentation for each method for more details.

## Output

Results of computations are returned to the REPL terminal by default. This can be controlled with the keywords `printstyle` and `output`.

- `output` contains the name of the file where results will be written. The default is `fermi.out`

- `printstyle` accepts four options

|    |    |
|----|----|
| **repl** | Print results to the REPL terminal |
| **file** | Print results to the file specified with the keyword `output` |
| **both** | Print results to the REPL and write it to the file specified with the keyword `output` |
| **none** | Does not print any results |

## Argument passing

In some situations, we may want to pass argument to the energy computation. There are two ways to do that in Fermi. 
1. Call the function directly, without using macros.
2. Use the `=>` or `<=` syntax within the `@energy` macro.

The code shown below computes a potential energy curve for the Helium dimer. First using functions directly
```julia
using Fermi
Rvals = [1.0 + 0.1*i for i = 0:10]
E = []
for r in Rvals
    mol = Molecule(molstring = """
        He 0.0 0.0 0.0
        He $r  0.0 0.0""")
    wfn = Fermi.MollerPlesset.RMP2(mol)
    push!(E, wfn.energy)
end
```
Now using the `@energy` macro
```julia
using Fermi
Rvals = [1.0 + 0.1*i for i = 0:10]
E = []
for r in Rvals
    mol = Molecule(molstring = """
        He 0.0 0.0 0.0
        He $r  0.0 0.0""")
    wfn = @energy mol => mp2
    push!(E, wfn.energy)
end
```
`@energy x,y,z,... => mp2` is equivalent to `Fermi.MollerPlesset.RMP2(x,y,z,...)`

## Interactive Usage

Interactive usage may be your best option for quick tasks, debugging, or if you are simply browsing the code. Moreover, this feature allows for usage within notebook environments such as [Jupyter](https://jupyter.org/try) and [Pluto](https://github.com/fonsp/Pluto.jl). Some objects, such as `Molecule`, `BasisSet` or `RHF`, can be printed directly on the terminal for some overview of their content.
```julia
julia> using Fermi
julia> using Fermi.Integrals

julia> mol = Molecule()
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 0   Multiplicity: 1   
Nuclear repulsion:    8.8880641737

julia> BasisSet("sto-3g", mol)
sto-3g Basis Set
Number of shells: 5
Number of basis:  7

O: 1s 2s 1p 
H: 1s 
H: 1s

julia> @set printstyle none;
julia> wfn = @energy rhf
 ⇒ Fermi Restricted Hartree--Fock Wave function
 ⋅ Basis:                  sto-3g
 ⋅ Energy:                 -74.965002894685
 ⋅ Occ. Spatial Orbitals:  5
 ⋅ Vir. Spatial Orbitals:  2
Convergence: ΔE => 0.00e+00 Dᵣₘₛ => 2.00e-10
```
While using Fermi interactively, one can use the macro `@get` to check the current value of an option keyword
```julia
julia> @set basis 6-31g
julia> @get basis
"6-31g"
```
or it can be used without arguments to check all user defined keywords
```julia
julia> @set {
  basis 6-31g
  df false
  diis false
}
julia> @get
┌─────────┬───────────────┐
│ Keyword │ Current Value │
├─────────┼───────────────┤
│    diis │         false │
│      df │         false │
│   basis │         6-31g │
└─────────┴───────────────┘
```
The macro `@reset` erases all options set
```julia
julia> @reset
julia> @get
No user defined keywords found.
```
Finally, the macro `@lookup` offers a quick way to look up keywords
```julia
julia> @lookup scf
┌──────────────┬───────────────┐
│      Keyword │ Current Value │
├──────────────┼───────────────┤
│  scf_max_rms │   1.00000e-09 │
│ scf_max_iter │            50 │
│    scf_guess │           gwh │
│      scf_alg │             1 │
│   scf_e_conv │   1.00000e-10 │
└──────────────┴───────────────┘
```

# Available methods
The following methods are currently implemented in Fermi

| Method    | Conventional | Density-Fitted | Single Precision |
|-----------|:------------:|:--------------:|:----------------:|
| RHF       |   ✔️          |    ✔️           |          ✖️       |
| UHF       |   ✔️          |    ✔️           |          ✖️       |
| RMP2      |   ✔️          |    ✔️           |    ✔️             |
| RCCSD     |   ✔️          |    ✔️           |    ✔️             |
| RCCSD(T)  |   ✔️          |    ✔️           |    ✔️             |

Only restricted reference methods are currently supported for correlated methods. All methods can use density fitting by setting `@set df true`. Moreover, JKFIT and RIFIT basis can be specified as
```julia
@set {
    jkfit cc-pvqz-jkfit
    rifit cc-pvqz-rifit
}
```
Single precision calculations are also possible using `@set precision single`.

# About
Fermi.jl is developed at the Center for Computational Quantum Chemistry at the University of Georgia under the supervision 
of [Dr. Justin M. Turney](https://github.com/jturney) and Prof. Henry F. Schaefer. For any questions, suggestions or if you want to participate in this project, please email [Gustavo Aroeira](https://github.com/gustavojra) (aroeira at uga.edu).