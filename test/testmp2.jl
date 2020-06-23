using Fermi
E = Fermi.Input.run("mp2_sto3g.dat")
@test isapprox(E,-0.04191336790304866;rtol=1E-10)
