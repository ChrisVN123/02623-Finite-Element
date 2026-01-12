using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")

function assembly_matrix(VX, VY, EToV)
    lam1 = 1 
    lam2 = 1

    N = size(EToV)[1]
    
    A = Matrix{Float64}(undef, N, N)
    b = Matrix{Float64}(undef, N, 1)
    
    for n in 1:N
        i, j, k = EToV[n, :]
        x_i, y_i = VX[i], VY[i]
        x_j, y_j = VX[j], VY[j]
        x_k, y_k = VX[k], VY[k]
        
        for r in 1:3
            Δ, abc = basfun(n, VX, VY, EToV)
            q̃ = 0.5 # FIX 
            q_r = abs(Δ)/3*q
            
            ĩ = 0 
            if r == 1 
                ĩ = i 
            elseif r == 2 
                ĩ = j 
            elseif r == 3 
                ĩ = k 
            end 
            
            b[ĩ] += q_r
            for s = 1:3
                j̃ = 0 
                if r == 1 
                    j̃ = i 
                elseif r == 2 
                    j̃ = j 
                elseif r == 3 
                    j̃ = k 
                end 

                k_rs = 1 / (4 abs(Δ)) * (lam1 * b_r * b_s + lam2 * c_r * c_s)
            end 
        end 


        q(n) = nothing
        b[i] := b[i] + q(n)
    
            for s in 1:3
                #Compute k(n)r,s
                if j ≥ i
                    #Look up the global numbers (i,j) and (x,y)-coordinates of the nodes in en.
                    a[i,j] := a[i,j] + k(n)r,s
                else
                    a[j,i] := a[j,i] + k(n)
                end
            end
    end