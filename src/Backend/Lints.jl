# Interface with Lints

using Lints
using LinearAlgebra

"""
    ao_kinetic(molecule::Fermi.Molecule, basis::String)

Computes AO basis kinetic energy ⟨μ|T̂|ν⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["T"]

where `helper` is bound to the desired molecule and basis set.
"""
function ao_kinetic(molecule::Molecule, basis::String; normalize=false)

    mol = mol_to_lints_molecule(molecule)
    bas = Lints.BasisSet(basis, mol)
    @lints begin
        T = Lints.make_T(bas; normalize = normalize)
    end
    return T
end

"""
    ao_overlap(molecule::Fermi.Molecule, basis::String)

Computes AO basis overlap ⟨p|q⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["S"]

where `helper` is bound to the desired molecule and basis set.
"""
function ao_overlap(molecule::Molecule, basis::String; normalize=false)

    mol = mol_to_lints_molecule(molecule)
    bas = Lints.BasisSet(basis, mol)
    @lints begin
        S = Lints.make_S(bas; normalize = normalize)
    end
    return Array(Hermitian(S))
end

"""
    ao_nuclear(molecule::Fermi.Molecule, basis::String)

Computes AO basis nuclear attraction ⟨μ|V̂|ν⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["V"]

where `helper` is bound to the desired molecule and basis set.
"""
function ao_nuclear(molecule::Molecule, basis::String; normalize=false)

    mol = mol_to_lints_molecule(molecule)
    bas = Lints.BasisSet(basis, mol)
    @lints begin
        V = Lints.make_V(bas; normalize=normalize)
    end
    return V
end

"""
    ao_eri(molecule::Fermi.Molecule, basis::String)

Computes AO basis electron repulsion integrals ⟨μν|Ô₂|ρσ⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["μ"]

where `helper` is bound to the desired molecule and basis set.
"""
function ao_eri(molecule::Molecule, basis::String; normalize=false)

    mol = mol_to_lints_molecule(molecule)
    bas = Lints.BasisSet(basis,mol)
    @lints begin
        I = Lints.make_ERI4(bas; normalize=normalize)
    end
    return I
end

"""
    df_ao_eri(molecule::Fermi.Molecule, basis::String)

Computes AO basis density fitted electron repulsion integrals ⟨μν|Ô₂|P⟩J(P,Q)^-1/2 integrals for the given basis and molecule.
Note that the returned integrals DO NOT need to be combined with the Coulomb metric J(P,Q). In common notation, this is B(Q,μ,ν).
Can be accessed at a higher level by calling
    
    helper["B"]

where `helper` is bound to the desired molecule and basis set.
"""
function df_ao_eri(molecule::Molecule, basis::String, aux::String; normalize=false)

    mol = mol_to_lints_molecule(molecule)
    bas = Lints.BasisSet(basis,mol)
    dfbas = Lints.BasisSet(aux,mol)
    @lints begin
        Pqp = Lints.make_ERI3(bas,dfbas; normalize=normalize)
        J = Lints.make_ERI2(dfbas; normalize=normalize)
    end
    Jh = Array(Hermitian(J)^(-1/2)) #sometimes Jh becomes complex slightly if J is not ~~exactly~~ hermitian 💔
    sz = Lints.nao(bas)
    for p=1:sz
        for q=1:sz
            auxP = Pqp[:,p,q]
            auxQ = Jh*auxP
            Pqp[:,p,q] .= auxQ
        end
    end
    return Pqp
end

function mol_to_lints_molecule(M::Molecule)
    natoms = length(M.atoms)
    zs = zeros(Int64,natoms)
    pos = fill(Float64[],natoms)
    for i=1:natoms
        zs[i] = M.atoms[i].Z
        pos[i] = collect(M.atoms[i].xyz)
    end
    Lints.Molecule(zs,pos)
end
