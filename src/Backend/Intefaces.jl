using Fermi.Geometry: Molecule
export @avogadro

macro avogadro(M)
    println("aa")
    return :(esc(Fermi.Geometry.get_string($M)))
end
