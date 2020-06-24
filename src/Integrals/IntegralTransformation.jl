"""
    Fermi.IntegralTransformation

Module to handle integral transformations from AO to MO.

**Functions:**

    get_eri   From a Wavefunction object, return a specified ERI array.
    get_fock  From a Wavefunction object, return the Fock matrix.

"""
module IntegralTransformation
using TensorOperations
using Fermi.Wavefunction

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
function get_eri(wfn::Wfn, eri_string::String; notation::String = "phys", fcn::Int = 0)

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
            o = 1+fcn:size(wfn.Cbo)[2]
            push!(C, wfn.Cbo[:,o])
            push!(S, wfn.nbeta)
        elseif s == 'O'
            o = 1+fcn:size(wfn.Cao)[2]
            push!(C, wfn.Cao[:,o])
            push!(S, wfn.nalpha)
        elseif s == 'v'
            push!(C, wfn.Cbv)
            push!(S, wfn.nvirb)
        elseif s == 'V'
            push!(C, wfn.Cav)
            push!(S, wfn.nvira)
        end
    end

    C1, C2, C3, C4 = C

    #gao = TBLIS.TTensor{eltype(wfn.ao_eri)}(wfn.ao_eri)
    #V = zeros(S[1],wfn.nmo,wfn.nmo,wfn.nmo)
    #v = TBLIS.TTensor{eltype(wfn.ao_eri)}(V)
    #c1 = TBLIS.TTensor{eltype(wfn.ao_eri)}(C1)
    #TBLIS.mul!(v,c1,gao,"ui","uvls","ivls")

    #V2 = zeros(S[1],S[2],wfn.nmo,wfn.nmo)
    #v2 = TBLIS.TTensor{eltype(wfn.ao_eri)}(V2)
    #c2 = TBLIS.TTensor{eltype(wfn.ao_eri)}(C2)
    #TBLIS.mul!(v2,c2,v,"va","ivls","ials")

    #V = zeros(S[1],S[2],S[3],wfn.nmo)
    #v = TBLIS.TTensor{eltype(wfn.ao_eri)}(V)
    #c3 = TBLIS.TTensor{eltype(wfn.ao_eri)}(C3)
    #TBLIS.mul!(v,c3,v2,"lj","ials","iajs")

    #V2 = zeros(S[1],S[2],S[3],S[4])
    #v2 = TBLIS.TTensor{eltype(wfn.ao_eri)}(V2)
    #c4 = TBLIS.TTensor{eltype(wfn.ao_eri)}(C4)
    #TBLIS.mul!(v2,c4,v,"sb","iajs","iajb")
    #Vnew = V2
    gao = wfn.ao_eri
    @tensoropt Vnew[i,a,j,b] := C4[σ,b]*C3[λ,j]*C2[ν,a]*C1[μ,i]*gao[μ,ν,λ,σ]

    if notation == "phys"
        @tensor Vnew[i,j,a,b] := Vnew[i,a,j,b]
    end

    return Vnew
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
function get_fock(wfn::Wfn; spin = "alpha")

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
