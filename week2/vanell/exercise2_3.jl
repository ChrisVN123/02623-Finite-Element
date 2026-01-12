using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")

function assembly(VX, VY, EToV, lam1, lam2, qt)
    """
    Algorithm 4 
    
    SHOULD BE SPARSE 
    Smart way? 
    """
    
    N = size(EToV)[1]
    M = length(VX)
    
    A = zeros(Float64, M, M)
    b = zeros(Float64, M)  # vector is simpler/safer

    for n ∈ 1:N
        i, j, k = EToV[n, :]
        Δ, abc = basfun(n, VX, VY, EToV)        
        
        for r ∈ 1:3
            q̃ = (qt[i] + qt[j] + qt[k]) / 3 # Can be moved out
            q_r_n = abs(Δ)/3 * q̃
            ĩ = EToV[n, r]
            
            b[ĩ] += q_r_n
            for s ∈ 1:3
                j̃ = EToV[n, s]

                b_r, b_s = abc[r, 2], abc[s, 2]
                c_r, c_s = abc[r, 3], abc[s, 3]
                
                k_rs = 1 / (4 * abs(Δ)) * (lam1 * b_r * b_s + lam2 * c_r * c_s)
                A[ĩ, j̃] += k_rs
            end 
        end 
    end

    A = sparse(A)
    return A, b 
end 