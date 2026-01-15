using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials
include("exercise2_1.jl")

function basfun(n, VX, VY, EToV)
    abc = Matrix{Float64}(undef, 3, 3)
    
    i, j, k = EToV[n, :]
    x1, y1 = VX[i], VY[i]
    x2, y2 = VX[j], VY[j]
    x3, y3 = VX[k], VY[k]

    Δ = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

    # (i,j,k) = (1, 2, 3)
    abc[1, 1] = x2*y3 - x3*y2   # ã₁
    abc[1, 2] = y2 - y3         # b̃₁
    abc[1, 3] = x3 - x2         # c̃₁

    # (i,j,k) = (2, 3, 1)
    abc[2, 1] = x3*y1 - x1*y3   # ã₂
    abc[2, 2] = y3 - y1         # b̃₂
    abc[2, 3] = x1 - x3         # c̃₂

    # (i,j,k) = (3, 1, 2)
    abc[3, 1] = x1*y2 - x2*y1   # ã₃ 
    abc[3, 2] = y1 - y2         # b̃₃
    abc[3, 3] = x2 - x1         # c̃₃
    return Δ, abc
end

function outernormal(n, k, VX, VY, EToV)
    i, j, h = EToV[n, :]
    x1, x2 = nothing, nothing
    y1, y2 = nothing, nothing
    
    if k == 1 
        x1, y1 = VX[i], VY[i]
        x2, y2 = VX[j], VY[j]
    elseif k == 2 
        x1, y1 = VX[j], VY[j]
        x2, y2 = VX[h], VY[h]
    elseif k == 3 
        x1, y1 = VX[h], VY[h] 
        x2, y2 = VX[i], VY[i]
    end 
    
    t1 = x2 - x1
    t2 = y2 - y1
    normal = 1 / (t1^2 + t2^2)^0.5
    
    return [t2, -t1] * normal
end

