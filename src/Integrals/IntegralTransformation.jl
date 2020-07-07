"""
    Fermi.IntegralTransformation

Module to handle integral transformations from AO to MO.

**Functions:**

    get_eri   From a Wavefunction object, return a specified ERI array.
    get_fock  From a Wavefunction object, return the Fock matrix.

"""
module IntegralTransformation
using TensorOperations
using Fermi

export get_eri
export get_fock

"""
    Fermi.IntegralTransformation.get_eri(wfn::Wfn, eri_string::String; notation::String = "phys")

From a Wavefunction object, return a specified ERI array.

**Arguments**

    wfn         Wavefunction object.
    eri_string  String with length 4 identifying the type of ERI. 
                Characters must be (o, O, v, V). Each indicating 
                Occupied and Virtual for ALPHA and beta.

**Kwargs**

    notation    {"chem", "phys"}
                Return the array in Chemist's or Physicists' notation.
                Default: "phys".

"""
function get_eri(wfn::wfT, eri_string::String; notation::String = "phys", fcn::Int = 0) where wfT <: Fermi.AbstractReferenceWavefunction

    # Size of eri_string must be 4.
    if sizeof(eri_string) != 4
        error("Invalid string given to Fermi.IntegralTransformation.get_eri: $eri_string")
    end

    # The computations below assumed Chemist's notation. Thus, the string is modified if Phys is requested
    if notation == "phys"
        eri_string = eri_string[[1,3,2,4]]
    end
    C = []
    S = []
    # Get C1, C2, C3, C4 for the integral transformation
    for s in eri_string
        if s == 'o'
            o = 1+fcn:size(Fermi.HartreeFock.RHF.Cbo(wfn))[2]
            push!(C, wfn.Cbo[:,o])
            push!(S, wfn.noccb)
        elseif s == 'O'
            o = 1+fcn:size(Fermi.HartreeFock.RHF.Cao(wfn))[2]
            push!(C, Fermi.HartreeFock.RHF.Cao(wfn)[:,o])
            push!(S, wfn.nocca)
        elseif s == 'v'
            push!(C, Fermi.Cbv(wfn))
            push!(S, wfn.nvirb)
        elseif s == 'V'
            push!(C, Fermi.HartreeFock.RHF.Cav(wfn))
            push!(S, wfn.nvira)
        end
    end

    C1, C2, C3, C4 = C

    gao = wfn.ERI
    T = typeof(gao)
    ## TODO: make intermediate tensors of same tensor type as gao
    nmo = Fermi.HartreeFock.RHF.nmo(wfn)
    Q1 = zeros(S[1],nmo,nmo,nmo)
    Fermi.contract!(Q1,C1,gao,"ivls","ui","uvls")
    Q2 = zeros(S[1],S[2],nmo,nmo)
    Fermi.contract!(Q2,C2,Q1,"ials","va","ivls")
    Q1 = nothing
    Q3 = zeros(S[1],S[2],S[3],nmo)
    Fermi.contract!(Q3,C3,Q2,"iajs","lj","ials")
    Q2 = nothing
    Q4 = zeros(S[1],S[2],S[3],S[4])
    Fermi.contract!(Q4,C4,Q3,"iajb","sb","iajs")

    if notation == "phys"
        Q4 = permutedims(Q4,(1,3,2,4))
    end
    return T(Q4)
end

"""
    Fermi.IntegralTransformation.get_fock(wfn::Wfn; spin::String = "alpha")

From a Wavefunction object, return the Fock matrix.

**Arguments**

    wfn  Wavefunction object.

**Kwargs**

    spin    {"alpha", "a", "up", "beta", "b", "down"}
            String indicating the spin of the Fock matrix. 
            Case insensitive. 

"""
function get_fock(wfn::wfT; spin = "alpha") where wfT <: Fermi.AbstractReferenceWavefunction

    if lowercase(spin) in ["alpha", "up", "a"]
        C  = wfn.Ca
        Co = wfn.Cao
    elseif lowercase(spin) in ["beta", "down", "b"]
        C  = wfn.Cb
        Co = wfn.Cbo
    else
        error("Invalid Spin option given to Fermi.IntegralTransformation.get_fock: $spin")
    end

    gao = wfn.ao_eri
    hao = wfn.hao

    @tensoropt (p=>100x, q=>100x, k=>x, μ=>100x, ν=>100x, λ=>100x, σ=>100x) begin
        f[p,q] := C[μ,p]*C[ν,q]*hao[μ,ν] 
        f[p,q] += 2*C[μ,p]*C[ν,q]*Co[λ,k]*Co[σ,k]*gao[μ,ν,λ,σ]
        f[p,q] -= C[μ,p]*C[ν,q]*Co[λ,k]*Co[σ,k]*gao[μ,λ,ν,σ]
    end

    return f
end


end # Module
