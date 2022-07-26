using GaussianBasis
using Molecules
using TensorOperations

function RMP2grad(x...)
    RMP2grad(Molecule(), x...)
end

function RMP2grad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(FermiException("Invalid or unsupported derivative type for RMP2: \"$dtype\""))
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.MollerPlesset.RMP2, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###