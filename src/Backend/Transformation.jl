"""
    Fermi.Transformation

this module manages all _complete_ integral transformations,
either disk-cached or in-core.
"""

module Transformation
using Fermi.DiskTensors
using TensorOperations
using LinearAlgebra

export tei_transform

function tei_transform(gao::Array{Float64,4},
                       C::Array{Float64,2},
                       name::String="default"
                      )
    tei_transform(gao,C,C,C,C,name)
end
function tei_transform(gao::DiskFourTensor,
                       C::Array{Float64,2},
                       name::String="default"
                      )
    tei_transform(gao,C,C,C,C,name)
end

"""
    tei_transform(
                  gao::Union{Array{Float64,4},Array{Float32,4}},
                  C1::Union{Array{Float64,2},Array{Float32,2}},
                  C2::Union{Array{Float64,2},Array{Float32,2}},
                  C3::Union{Array{Float64,2},Array{Float32,2}},
                  C4::Union{Array{Float64,2},Array{Float32,2}},
                  name::String
                  )
Transform `gao` to the MO basis defined by C1,C2,C3,C4. 
"""
function tei_transform(
                       gao::Union{Array{Float64,4},Array{Float32,4}},
                       C1::Union{Array{Float64,2},Array{Float32,2}},
                       C2::Union{Array{Float64,2},Array{Float32,2}},
                       C3::Union{Array{Float64,2},Array{Float32,2}},
                       C4::Union{Array{Float64,2},Array{Float32,2}},
    name::String,
)
    """
    general spin orbital transformation
     
    """
    T = typeof(gao)
    norb = size(C1)[1]
    d = size(C1)[1]
    d1 = size(C1)[2]
    d2 = size(C2)[2]
    d3 = size(C3)[2]
    d4 = size(C4)[2]
    rr = UnitRange(1, norb)
    rr1 = UnitRange(1, d1)
    rr2 = UnitRange(1, d2)
    rr3 = UnitRange(1, d3)
    rr4 = UnitRange(1, d4)
    R = collect(rr)
    R1 = collect(rr1)
    R2 = collect(rr2)
    R3 = collect(rr3)
    R4 = collect(rr4)
    temp = zeros(d, d, d, d4)
    #Quarter transform 1
    #(μ,ν,λ,σ) -> (μ,ν,λ,b)
    @tensoropt begin
        temp[μ,ν,λ,b] = C4[σ,b]*gao[μ,ν,λ,σ]
    end
    #Quarter transform 2
    #(μ,ν,λ,b) -> (μ,ν,j,b)
    temp2 = zeros(d, d, d3, d4)
    @tensoropt begin
        temp2[μ,ν,j,b] = C3[λ,j]*temp[μ,ν,λ,b]
    end
    #Quarter transform 3
    #(μ,ν,j,b) -> (μ,a,j,b)
    temp = zeros(d, d2, d3, d4)
    @tensoropt begin
        temp[μ,a,j,b] = C2[ν,a]*temp2[μ,ν,j,b]
    end
    #Quarter transform 4
    #(μ,a,j,b) -> (i,a,j,b)
    temp2 = zeros(d1, d2, d3, d4)
    @tensor begin
        temp2[i,a,j,b] = C1[μ,i]*temp[μ,a,j,b]
    end
    return temp2
end
function tei_transform(
    gao::DiskFourTensor,
    C1::Union{Array{Float64,2},Array{Float32,2}},
    C2::Union{Array{Float64,2},Array{Float32,2}},
    C3::Union{Array{Float64,2},Array{Float32,2}},
    C4::Union{Array{Float64,2},Array{Float32,2}},
    name::String,
)
    """
    general spin orbital transformation
     
    """
    T = eltype(gao)
    norb = size(C1)[1]
    d = size(C1)[1]
    d1 = size(C1)[2]
    d2 = size(C2)[2]
    d3 = size(C3)[2]
    d4 = size(C4)[2]
    rr = UnitRange(1, norb)
    rr1 = UnitRange(1, d1)
    rr2 = UnitRange(1, d2)
    rr3 = UnitRange(1, d3)
    rr4 = UnitRange(1, d4)
    R = collect(rr)
    R1 = collect(rr1)
    R2 = collect(rr2)
    R3 = collect(rr3)
    R4 = collect(rr4)
    temp = DiskFourTensor("/tmp/jues.$name.temp.0", T, d, d, d, d4, "w")
    blockfill!(temp, 0.0)
    #Quarter transform 1
    #(μ,ν,λ,σ) -> (μ,ν,λ,b)
    ocache = zeros(d,d,d4)
    icache = zeros(size(gao[1,:,:,:]))
    for μ in R
        #for ν in R
        ocache .= 0.0#zeros(d,d4)
        icache .= gao[μ,:,:,:]
        @tensoropt begin
            ocache[ν,λ,b] = C4[σ,b]*icache[ν,λ,σ]
        end
        temp[μ,:,:,:] .= ocache[:,:,:]
        #end
    end
    #Quarter transform 2
    #(μ,ν,λ,b) -> (μ,ν,j,b)
    temp2 = DiskFourTensor("/tmp/jues.$name.temp2.0", T, d, d, d3, d4, "w")
    blockfill!(temp2, 0.0)
    for μ in R
        for ν in R
            ocache = zeros(d3,d4)
            icache = temp[μ,ν,:,:]
            @tensoropt begin
                ocache[j,b] = C3[λ,j]*icache[λ,b]
            end
            temp2[μ,ν,:,:] = ocache[:,:]
        end
    end
    #Quarter transform 3
    #(μ,ν,j,b) -> (μ,a,j,b)
    temp = DiskFourTensor("/tmp/jues.$name.temp.0", T, d, d2, d3, d4, "w")
    blockfill!(temp, 0.0)
    for b in R4
        for j in R3
            ocache = zeros(d, d2)
            icache = temp2[:, :, j, b]
            for a in R2
                for μ in R
                    for ν in R
                        ocache[μ, a] += C2[ν, a] * icache[μ, ν]
                    end
                end
            end
            temp[:, :, j, b] = ocache[:, :]
        end
    end
    #Quarter transform 4
    #(μ,a,j,b) -> (i,a,j,b)

    temp2 = DiskFourTensor("/tmp/jues.$name.temp2.0", T, d1, d2, d3, d4, "w")
    blockfill!(temp2, 0.0)
    for b in R4
        for j in R3
            icache = temp[:, :, j, b]
            ocache = zeros(d1, d2)
            for a in R2
                for i in R1
                    for μ in R
                        #temp2[i,a,j,b] += C1[μ,i]*temp[μ,a,j,b]
                        ocache[i, a] += C1[μ, i] * icache[μ, a]
                    end
                end
            end
            temp2[:, :, j, b] = ocache
        end
    end
    return temp2
end

end # Module
