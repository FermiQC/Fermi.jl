# Fermi.jl
A Julia module for *ab initio* quantum chemistry computations.

# Installation

Fermi.jl is a registred Julia package and, as such, it can be obtained using the standard Julia package manager. From the Julia terminal, use the `]` to move to the pkg manager
<pre><code><b style="color: LimeGreen">julia></b> # This is the standard Julia terminal, hit ]
<b style="color: DodgerBlue">(@v1.6) pkg> </b> # This is the package manager! Hit back space to leave this mode</code></pre>
Next add Fermi to the current enviroment.
<pre><code><b style="color: DodgerBlue">(@1.6) pkg></b> Add Fermi</pre></code>
All the dependencies are going to be downloaded and installed and the code should be ready to work. To test the package you can run
<pre><code><b style="color: DodgerBlue">(@1.6) pkg></b> test Fermi</pre></code>

## Trouble-shooting
The code is built and tested for the latest Julia version on Ubuntu and MacOS. The most fragile step that can lead to problems while building Fermi is the integral library `libcint`. Namely, `libcint` must be able to find `BLAS` in your computer. For Linux, `BLAS` can be easily installed as
```
sudo apt-get install libblas-dev liblapack-dev
```
You might want to check the `build.jl` file that contains the instructions to fetch and install `libcint`. This file is in the source code of Fermi. Julia store source codes in the `.julia` folder. If you used `add Fermi` to install the package, the source code should be located at `.julia/packages/Fermi`. Please also refer to the `libcint` [github page](https://github.com/sunqm/libcint) for more details on the dependecies there. 

# Usage

Fermi.jl can be used interactively through the Julia terminal, or you can write a Julia script which will act as the traditional input file present in other quantum chemistry packages. For example, a minimal script to run a RHF computation on a water molecule can be written as
```
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

Interactevely, the functions calls are all the same. First you load Fermi in
<pre><code><b style="color: LimeGreen">julia></b> using Fermi</code></pre>
Next you set the molecule and the options desired.
<pre><code><b style="color: LimeGreen">julia></b> @molecule {
  O        1.2091536548      1.7664118189     -0.0171613972
  H        2.1984800075      1.7977100627      0.0121161719
  H        0.9197881882      2.4580185570      0.6297938830
}
<b style="color: LimeGreen">julia></b> @set {
    basis cc-pvdz
    df true
}</pre></code>
Finally, run the desired energy computations
<pre><code><b style="color: LimeGreen">julia></b> @energy mp2</code></pre>

# Available methods
The following methods are currently implemented in Fermi
<table align="center">
  <tr>
    <th>Method</th>
    <th>Conventional</th>
    <th>Density-Fitted</th>
    <th>Single Precision</th>
  </tr>

  <tr align="center">
    <td>RHF</td>
    <td>&#10004</td>
    <td>&#10004</td>
    <td>&#10060</td>
  </tr>

  <tr align="center">
    <td>RMP2</td>
    <td>&#10004</td>
    <td>&#10004</td>
    <td>&#10004</td>
  </tr>

  <tr align="center">
    <td>RCCSD</td>
    <td>&#10004</td>
    <td>&#10004</td>
    <td>&#10004</td>
  </tr>

  <tr align="center">
    <td>RCCSD(T)</td>
    <td>&#10004</td>
    <td>&#10004</td>
    <td>&#10004</td>
  </tr>
</table>

Only restricted reference methods are currently supported. However, all methods can be run using density fitting by setting `@set df true`. Moreover, JKFIT and RIFIT basis can be specified as
```
@set {
    jkfit cc-pvqz-jkfit
    rifit cc-pvqz-rifit
}
```
Single precision calculations are also possible using `@set precision single`.

# About
Fermi.jl is developed at the Center for Computational Quantum Chemistry at the University of Georgia under the supervision 
of [Dr. Justin M. Turney](https://github.com/jturney) and Prof. Henry F. Schaefer. For any questions, suggestions or if you want to participate in this project, please email [Gustavo Aroeira](https://github.com/gustavojra) (aroeira at uga.edu).