### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 4e019dd0-9a70-4034-915b-baf59e3416e5
using Fermi

# ╔═╡ e45bf90a-9cfe-46e8-bb0f-1c42723f1389
using TensorOperations

# ╔═╡ 819f62ec-da7c-11eb-300f-333f486858d2
html"""
<p align="center">
  <img src="https://raw.githubusercontent.com/FermiQC/Fermi.jl/dae4255e264c48e98810b22da3e7ee118ec649cd/docs/src/assets/logo.svg" width="400" alt=""/>
</p>
<h1 align="center"> <i> JuliaCon 2021 </i></h1>
<p align="center" style="font-size:25px"> Gustavo Aroeira <br/></p>
<p align="center"><b>Center for Computational Quantum Chemistry</b><br/>
University of Georgia, Athens, GA<br/>
United States of America
</p>
"""

# ╔═╡ e1bd2ce0-a958-4460-b48f-620db0743e45
md"""
## 1. Quantum Chemistry
The main goal in quantum chemistry is to solve the electronic Schrödinger equation for many-electron systems

$\hat{H}|\Psi\rangle = E|\Psi\rangle$

A first approximation to the fermionic wave function, known as *mean-field*, is represented as a Slater determinant

$|\Phi_0\rangle = \mathcal{A}[\phi_1(\mathbf{x_1})\phi_2(\mathbf{x_2})\phi_3(\mathbf{x_3})...\phi_N(\mathbf{x_N})]$

where $\mathcal{A}$ is the antisymmetrizer operator and $\phi_i$ are spin-orbitals, that is, one-electron wave functions used to build the many-electron *ansatz*. Using the *basis set approximation*, we represent the orbitals as an expansion of $M$ pre-selected basis functions

$\phi_i = \sum_\mu^M C_{\mu i} \chi_\mu$

This is known as the **L**inear **C**ombination of **A**tomic **O**rbitals (LCAO) as the basis functions are designed to resemble atomic orbitals. Finally, more precise, i.e. correlated, wave functions are built as combinations of different Slater determinants (they differ by the orbitals within each one).

$|\Psi\rangle = \sum_i c_i|\Phi_i\rangle$

or, in cluster notation

$|\Psi\rangle =  e^{T}|\Phi_0\rangle$

**Fermi** is a quantum chemistry package in Julia specialized in wave function methods. Currently we have implemented the following methods

- Restricted Hartree-Fock
- Møller-Plesset Perturbation Theory
- Coupled Cluster Theory

Moreover, Fermi aims to serve as development platform for molecular quantum mechanics in Julia
"""

# ╔═╡ ff5a766b-716b-47dc-a94f-1f9c3f74ac8b
md"""
## 2. Molecular Integrals

Basic quantities required in all quantum chemistry computations

###### Overlap (S)

$S_{ij} = \int \phi_i^* \phi_j d\tau$

###### Kinetic (T)

$T_{ij} = \int \phi_i^* \left( -\frac{1}{2} \nabla^2 \right) \phi_j d\tau$

###### Nuclear (V)

$V_{ij} = \sum_k \int \phi_i^* \frac{Z_k}{r} \phi_j d\tau$

###### Eletron repulsion integral (ERI)

$G_{ijkl} = \int \phi_i^* \phi_j^* \frac{1}{r_{12}} \phi_k \phi_l d\tau$

Fermi utilizes the library [libcint](https://github.com/sunqm/libcint) to obtain molecular integrals.
"""

# ╔═╡ 67892086-7643-40a0-9e5e-eb56607412de
begin
@set {
	diis false
	oda false
}
nothing
end

# ╔═╡ e8585228-0be6-41dd-8434-045ed69d578c
md"""
The module `Fermi.Libcint` has minimal wraps around libcint functions
"""

# ╔═╡ 3a325265-7699-4f99-a850-207ec00f9c26
names(Fermi.Libcint)

# ╔═╡ 8dced24e-6cc6-46aa-bd50-658907190987
md"""
A more convinient way is to use the `IntegralHelper`.

First, define a molecule and basis set
"""

# ╔═╡ 01228770-9ad9-4634-9aba-afc798ae3e8c
@molecule {
	H 0.0 0.0 0.0
	H 0.76 0.0 0.0
}

# ╔═╡ 597aa086-dd0c-4ef9-b0ef-ab294b89b3a2
@set basis sto-3g

# ╔═╡ bd1dca4c-f126-4353-b09a-e2a641712a09
I = Fermi.Integrals.IntegralHelper()

# ╔═╡ c37cef53-ddfa-4bce-9328-21cde32bf8da
I["S"]

# ╔═╡ 73551bfb-c7c6-4f5b-a931-92c23b7feafd
md"""
## 3. Multiple dispatch in Fermi

`IntegralHelper` guides how amplitudes are updated in a CCSD computation
"""

# ╔═╡ b03d9e3f-8ef5-4e48-b2aa-f6b4295ac228
html"""
<p align="center">
  <img src="https://raw.githubusercontent.com/FermiQC/Fermi.jl/dd6fd813d4143565146a36c7b84a49058908c7f7/docs/src/assets/ccflow.svg" width="600" alt=""/>
"""

# ╔═╡ fb42172a-33df-4b03-88d5-6db7224f132b
md"""
## 4. Tensor Contractions

In methods such as Coupled Cluster, tensor contractions appears in terms like

$t_{ij}^{ab} \leftarrow -2t_{ik}^{ac}\cdot t_{lj}^{bd}\cdot v_{kc}^{ld}$

The package [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl) is used to handle expressions of this type.

`FermiMDArray` is a thin wrap around a Julia array that allow us to use [TBLIS](https://github.com/devinamatthews/tblis) as a tensor contraction tool.
"""

# ╔═╡ d5942379-2e27-43ca-920e-889a64d32269
begin
	o_size = 10
	v_size = 100
	T2 = FermiMDrand(o_size, o_size, v_size, v_size)
	newT2 = FermiMDrand(o_size, o_size, v_size, v_size)
	Vovov = FermiMDrand(o_size, v_size, o_size, v_size)
	size(Vovov)
end

# ╔═╡ d0c5639e-49f7-45c0-a19e-d1da88b81ce0
t = @elapsed begin
	@tensoropt newT2[i,j,a,b] = -2*T2[i,k,a,c]*T2[l,j,b,d]*Vovov[k,c,l,d]
end

# ╔═╡ e89556da-e726-47bf-ae91-4312b24095b1
@set tblis true

# ╔═╡ Cell order:
# ╟─819f62ec-da7c-11eb-300f-333f486858d2
# ╟─e1bd2ce0-a958-4460-b48f-620db0743e45
# ╟─ff5a766b-716b-47dc-a94f-1f9c3f74ac8b
# ╟─67892086-7643-40a0-9e5e-eb56607412de
# ╠═4e019dd0-9a70-4034-915b-baf59e3416e5
# ╟─e8585228-0be6-41dd-8434-045ed69d578c
# ╠═3a325265-7699-4f99-a850-207ec00f9c26
# ╟─8dced24e-6cc6-46aa-bd50-658907190987
# ╠═01228770-9ad9-4634-9aba-afc798ae3e8c
# ╠═597aa086-dd0c-4ef9-b0ef-ab294b89b3a2
# ╠═bd1dca4c-f126-4353-b09a-e2a641712a09
# ╠═c37cef53-ddfa-4bce-9328-21cde32bf8da
# ╟─73551bfb-c7c6-4f5b-a931-92c23b7feafd
# ╟─b03d9e3f-8ef5-4e48-b2aa-f6b4295ac228
# ╟─fb42172a-33df-4b03-88d5-6db7224f132b
# ╠═e45bf90a-9cfe-46e8-bb0f-1c42723f1389
# ╠═d5942379-2e27-43ca-920e-889a64d32269
# ╠═d0c5639e-49f7-45c0-a19e-d1da88b81ce0
# ╠═e89556da-e726-47bf-ae91-4312b24095b1
