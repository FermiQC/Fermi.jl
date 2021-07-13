### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 416d9c62-d2ef-11eb-3a40-450aee175208
using Fermi.GaussianBasis

# ╔═╡ 71915091-164d-440f-931e-7abaa1764e77
using Plots

# ╔═╡ 913f76eb-3087-4a85-ae69-5dd7565964e9
using Fermi.Geometry

# ╔═╡ aec67321-b943-478b-9765-8cb95c0f9d06
using Fermi.Integrals

# ╔═╡ f41f12be-3755-4f73-a183-497ced9afa24
using Fermi.HartreeFock

# ╔═╡ 42b34173-623d-4361-a2ba-b55b0cdfd69f
using Fermi

# ╔═╡ 1ad27f33-fe6c-498a-b4e8-89948dd1f94b
html"""
<p align="center">
  <img src="https://raw.githubusercontent.com/FermiQC/Fermi.jl/dae4255e264c48e98810b22da3e7ee118ec649cd/docs/src/assets/logo.svg" width="400" alt=""/>
</p>
<h1 align="center"> <i> Basic structures </i></h1>
<p align="center" style="font-size:25px"> Gustavo Aroeira <br/></p>
<p align="center"><b>Center for Computational Quantum Chemistry</b><br/>
University of Georgia, Athens, GA<br/>
United States of America
</p>
"""

# ╔═╡ c51d70f8-9065-41f9-aea3-73f7c2aba7e7
md"""
## 1. Goal
The main goal in quantum chemistry is to solve the electronic Schrödinger equation for many-electron systems

$\hat{H}|\Psi\rangle = E|\Psi\rangle$

A first approximation to the fermionic wave function is represented as a Slater determinant

$|\Psi_0\rangle = \mathcal{A}[\phi_1(\mathbf{x_1})\phi_2(\mathbf{x_2})\phi_3(\mathbf{x_3})...\phi_N(\mathbf{x_N})]$

where $\mathcal{A}$ is the antisymmetrizer operator and $\phi_i$ are spin-orbitals, that is, one-electron wave functions used to build the many-electron *ansatz*.

The first task is, therefore, to find optimized one-electron functions (orbitals). For atoms, one can do it numerically in a self-consistent fashion. For molecules however, this becomes untractable. Thus, we introduce the *basis set approximation*, where we represent the orbitals as an expansion of $M$ pre-selected basis functions

$\phi_i = \sum_\mu^M C_{\mu i} \chi_\mu$

This is known as the **L**inear **C**ombination of **A**tomic **O**rbitals (LCAO) as the basis functions are designed to resemble atomic orbitals.
"""

# ╔═╡ dc46fbbd-b96e-4fdd-9300-fa5ca87e86de
md"""
## 2. Basis Functions

Analytical solution to the ground state wave function of the hydrogen atom has the shape of a Slater orbital

$\Psi_\text{H} \sim e^{-r}$

Computing integrals and derivatives of slater type functions is too complicated. Thus, gaussian functions are employed instead. Each basis function in our basis set is a combination of gaussian type functions is to mimicking slater type functions.

$\chi = \sum_i N_i G_i(r,\theta,\phi)$
$G_i(r,\theta,\phi) = Y_{l,m_l}(\theta,\phi) r^l e^{-\zeta_i r^2}$

Each basis function $\chi$ is characterized by three quantities:

1. Angular momentum number $l$
2. A set of contraction coefficients $N_i$
3. A set of exponent values $\zeta_i$
"""

# ╔═╡ 97282505-47ee-4011-8198-bd4a7f8f75bf
md"""
##### Let's build a basis function using Fermi
Import the GaussianBasis module
"""

# ╔═╡ 2f8ca46a-3a92-4d0f-a28b-7e47fc577c80
begin
	l = 1
	exp_vals = [10.0, 5.0, 0.6]
	coef_vals = [0.2, 0.6, 0.7]
	BasisFunction(l, coef_vals, exp_vals)
end

# ╔═╡ 51444ecc-7031-4f9b-bb67-733bb6def9dc
@bind bname html"""Select basis set   <select>
    <optgroup label="Pople">
      <option value="sto-3g">STO-3g</option>
      <option value="sto-6g">STO-6g</option>
	  <option value="6-31g">6-31g</option>
    </optgroup>
    <optgroup label="Dunning">
      <option value="cc-pvdz">cc-pVDZ</option>
      <option value="cc-pvtz">cc-pVTZ</option>
      <option value="cc-pvqz">cc-pVQZ</option>
      <option value="cc-pv5z">cc-pV5Z</option>
    </optgroup>
  </select>"""

# ╔═╡ 3a9d4c13-9d77-4deb-add2-b54e29d59934
begin
	# Bohr radius
	a0 = 0.529177210903
	
	# Define Slater function
	slater(r) = exp(-r/a0)
	
	# Gaussian function
	function contracted_gaussian(r,e,c)
	out = zeros(length(r))
	for i = eachindex(e)
		ζ = e[i]
		C = c[i]
		out += C*exp.(-ζ.*r.^2 ./a0^2)
	end
	return out ./ maximum(out)
	end
	
	# Plot range
	rvals = collect(range(0,2,step=0.01))
	
	bfp, = GaussianBasis.read_basisset(bname, "H")
	plot(rvals, slater.(rvals), label="Slater orbital")
	plot!(rvals, contracted_gaussian(rvals, bfp.exp, bfp.coef), label=bname)
	xlabel!("Distance from the nucleus (Å)")
	ylabel!("ϕ(r)")
end

# ╔═╡ 20b7f076-7177-41fe-b0d6-abd1b9211360
md"""
## 3. Basis Set
"""

# ╔═╡ 57f7f7c7-8a84-4fd0-b3f9-f49d375bf37c
Atom("H", 1, (0.0,0.0,0.0))

# ╔═╡ 7388eaa1-011c-43df-9986-2f59ba259873
atoms = [Atom("H", 1, (0.0,0.0,0.0)),
	     Atom("H", 1, (0.76,0.0,0.0))]

# ╔═╡ 77186eea-5b85-40c9-91bf-f48707752bb7
mol = Molecule(atoms, 0, 1)

# ╔═╡ ee4dbc82-6258-413c-aca5-1bd1e6ee02fc
begin
	# Create basis functions
	s = BasisFunction(0, [0.5215367271], [0.122])
	p = BasisFunction(1, [1.9584045349], [0.727])
	
	# Map atoms to basis functions
	shells = Dict(
	atoms[1] => [s],
	atoms[2] => [s,p]
	)
	bset = BasisSet(mol, "Hs/Hsp", shells)
end

# ╔═╡ f17bc090-cad6-427a-ace9-7ed921f8368d
md"""
## 4. Integrals

In order to optimize orbitals several integrals over our basis functions are needed.
For example, the overlp integral

$S_{ij} = \int \phi_i^* \phi_j d\tau = C_{\mu i}^*C_{\nu j}\int \chi_\mu^* \chi_\nu d\tau$

Fermi uses the [libcint](https://github.com/sunqm/libcint) library to compute integrals.
"""

# ╔═╡ 15cfa11b-78e1-456c-82f2-02bcd1f8c6c8
I = IntegralHelper(bset)

# ╔═╡ 3685dd55-cb6a-4e0d-8a64-a76fd29d8e49
I["S"]

# ╔═╡ b6e3d8a9-a125-417c-8af1-145b9050d406
md"""
## 5. Hartree-Fock Orbitals
Having an integral object built, a Hartree--Fock computation can be run to obtain optimized orbitals
"""

# ╔═╡ 6fb41a7c-2ae6-4b98-ba5a-fe122512364a
wfn = RHF(I)

# ╔═╡ 90087552-51d4-42ef-8c3e-e5a61eec4fe1
wfn.orbitals.C

# ╔═╡ 683b958d-7959-484e-baf1-9803094dd994
md"""
## 6. Macros to make life easier

Fermi has a `Options` module that is referred to when some information is not specified. Moreover, users can use the macro `@molecule` and `@set` to change keywords and run computations using the desired settings.
"""

# ╔═╡ 165a0055-0195-42d8-8dd2-a3d07fe67f61
begin
	@molecule {
		H 0.0 0.0 0.0
		H 0.76 0.0 0.0
	}
	
	@set {
		basis 6-31g
		diis false
	}

	@energy rhf
end

# ╔═╡ 9ab18d9f-28a2-43ae-a9df-012223dea282
begin
	Ecc = []
	Emp2 = []
	Eref = []
	dist = collect(range(0.4,3.0,step=0.1))
	for d in dist
		Fermi.Options.set("molstring", 
			"""
			H 0.0 0.0 0.0
			H $d 0.0 0.0
			""")
		
		Iao = IntegralHelper()
		ref = @energy Iao=>rhf
		Imo = IntegralHelper(orbitals=ref.orbitals)
		cc_wfn = @energy Imo,Iao=>ccsd
		push!(Eref, ref.energy)
		push!(Emp2, cc_wfn.guessenergy)
		push!(Ecc, cc_wfn.energy)
	end
end

# ╔═╡ 8b86034b-c223-4c29-89ef-b03c24e736df
@bind rhf_plot html"""
<input type="checkbox" id="rhf">
<label for="rhf">RHF</label><br>
"""

# ╔═╡ 8f15233b-3724-437f-a78c-460f69a41da6
	@bind mp2_plot html"""
	<input type="checkbox" id="mp2">
	<label for="mp2">MP2</label><br>
	"""

# ╔═╡ 0d4ae0b1-6404-4471-818e-6cee43161b64
	@bind cc_plot html"""
	<input type="checkbox" id="mp2">
	<label for="cc">CCSD</label><br>
	"""

# ╔═╡ 850513be-8a24-4d4e-9d79-43062efaa50f
begin
	pt = plot(ylim=[-1.15, -0.80],legend=:topleft)
	xlabel!("H-H distance (Å)")
	ylabel!("Energy (Hartrees)")
	rhf_plot ? scatter!(dist, Eref, label="RHF", color=:steelblue2) : nothing
	mp2_plot ? scatter!(dist, Emp2, label="MP2",color=:seagreen) : nothing
	cc_plot ? scatter!(dist, Ecc, label="CCSD",color=:crimson) : nothing
	pt
end

# ╔═╡ Cell order:
# ╟─1ad27f33-fe6c-498a-b4e8-89948dd1f94b
# ╟─c51d70f8-9065-41f9-aea3-73f7c2aba7e7
# ╟─dc46fbbd-b96e-4fdd-9300-fa5ca87e86de
# ╟─97282505-47ee-4011-8198-bd4a7f8f75bf
# ╠═416d9c62-d2ef-11eb-3a40-450aee175208
# ╟─71915091-164d-440f-931e-7abaa1764e77
# ╠═2f8ca46a-3a92-4d0f-a28b-7e47fc577c80
# ╟─51444ecc-7031-4f9b-bb67-733bb6def9dc
# ╟─3a9d4c13-9d77-4deb-add2-b54e29d59934
# ╟─20b7f076-7177-41fe-b0d6-abd1b9211360
# ╠═913f76eb-3087-4a85-ae69-5dd7565964e9
# ╠═57f7f7c7-8a84-4fd0-b3f9-f49d375bf37c
# ╠═7388eaa1-011c-43df-9986-2f59ba259873
# ╠═77186eea-5b85-40c9-91bf-f48707752bb7
# ╠═ee4dbc82-6258-413c-aca5-1bd1e6ee02fc
# ╟─f17bc090-cad6-427a-ace9-7ed921f8368d
# ╠═aec67321-b943-478b-9765-8cb95c0f9d06
# ╠═15cfa11b-78e1-456c-82f2-02bcd1f8c6c8
# ╠═3685dd55-cb6a-4e0d-8a64-a76fd29d8e49
# ╟─b6e3d8a9-a125-417c-8af1-145b9050d406
# ╠═f41f12be-3755-4f73-a183-497ced9afa24
# ╠═6fb41a7c-2ae6-4b98-ba5a-fe122512364a
# ╠═90087552-51d4-42ef-8c3e-e5a61eec4fe1
# ╟─683b958d-7959-484e-baf1-9803094dd994
# ╠═42b34173-623d-4361-a2ba-b55b0cdfd69f
# ╠═165a0055-0195-42d8-8dd2-a3d07fe67f61
# ╟─9ab18d9f-28a2-43ae-a9df-012223dea282
# ╟─8b86034b-c223-4c29-89ef-b03c24e736df
# ╟─8f15233b-3724-437f-a78c-460f69a41da6
# ╟─0d4ae0b1-6404-4471-818e-6cee43161b64
# ╟─850513be-8a24-4d4e-9d79-43062efaa50f
