@set basis cc-pvdz
@test (@get basis) == "cc-pvdz"
@reset
@test (@get basis) == "sto-3g"
