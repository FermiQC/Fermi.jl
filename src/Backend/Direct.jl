module Direct
using Fermi.Wavefunction
using PyCall
using TensorOperations
export ao
export func_shell
export ao_to_mo_shell

function ao_function(wfn::DirectWfn, p, q, r, s)
    """
    Computes a single ERI in AO basis
    Very wasteful (throws out rest of shell)
    use for debugging
    """
    shell = ao_shell(wfn, p, q, r, s)
    return shell[po, qo, ro, so]
end
function func_shell(wfn::DirectWfn)
    """
    Matches basis function number and corresponding shell number
    usage: output[shellno] -> [fn1,fn2,...,fnn]
    """
    basis = wfn.basis
    nao = basis.nao()
    output = Dict()
    for i = 1:1:nao
        x = basis.function_to_shell(i - 1) + 1
        if haskey(output, x)
            push!(output[x], i)
        else
            output[x] = [i]
        end
    end
    println(output)
    return output
end
function ao_shell(wfn::DirectWfn, p, q, r, s)
    """
    Computes the AO eri shell (pq|rs)
    """
    basis = wfn.basis
    mints = wfn.mints
    ps = basis.function_to_shell(p - 1)
    qs = basis.function_to_shell(q - 1)
    rs = basis.function_to_shell(r - 1)
    ss = basis.function_to_shell(s - 1)

    po = p - basis.shell_to_ao_function(ps)# + 1
    qo = q - basis.shell_to_ao_function(qs)# + 1
    ro = r - basis.shell_to_ao_function(rs)# + 1
    so = s - basis.shell_to_ao_function(ss)# + 1
    shell = mints.ao_eri_shell(ps, qs, rs, ss).to_array()
    return shell
end
"""
    direct_ao_contract

RHF only
contraction to N^3 (iν|λσ) given specific i
"""
function direct_ao_contract(i,C::Array{Float64,2},mints::PyObject, basis::PyObject)
    nsh = basis.nshell()
    nao = basis.nbf()
    iνλσ = zeros(nao,nao,nao)
    Ci = C[:,i]
    for s in 1:nsh
        for r in 1:nsh
            for q in 1:nsh
                for p in 1:nsh
                    sf = basis.shell_to_basis_function(s-1) + 1
                    rf = basis.shell_to_basis_function(r-1) + 1
                    qf = basis.shell_to_basis_function(q-1) + 1
                    pf = basis.shell_to_basis_function(p-1) + 1
                    shell = mints.ao_eri_shell(p-1,q-1,r-1,s-1).np
                    pn, qn, rn, sn = size(shell)
                    pn -= 1
                    qn -= 1
                    rn -= 1
                    sn -= 1
                    #println("$pf $pn $qf $qn $rf $rn $sf $sn")
                    _Ci = Ci[pf:pf+pn]
                    @tensoropt begin
                        temp[ν,λ,σ] := _Ci[μ]*shell[μ,ν,λ,σ]
                    end
                    ##println("$p $q $r $s")
                    iνλσ[qf:qf+qn,rf:rf+rn,sf:sf+sn] += temp
                end
            end
        end
    end
    return iνλσ

end
end#module
