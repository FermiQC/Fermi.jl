using GaussianBasis
using Molecules
using TensorOperations

function UHFgrad(x...)
    UHFgrad(Molecule(), x...)
end

function UHFgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(FermiException("Invalid or unsupported derivative type for UHF: \"$dtype\""))
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.HartreeFock.UHF, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###