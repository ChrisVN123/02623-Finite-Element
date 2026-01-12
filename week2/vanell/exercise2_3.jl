using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")

println("EXERCISE 2.3")

function k_rs(r, s, lam1, lam2)
    
end 

function assembly_matrix(VX, VY, EToV, lam1, lam2, qt)
    N = size(EToV)[1]
    M = length(VX)
    
    A = zeros(Float64, M, M)
    b = zeros(Float64, M)  # vector is simpler/safer

    for n in 1:N
        i, j, k = EToV[n, :]
        Δ, abc = basfun(n, VX, VY, EToV)        
        
        for r in 1:3
            q̃ = (qt[i] + qt[j] + qt[k]) / 3 
            q_r_n = abs(Δ)/3 * q̃
            ĩ = EToV[n, r]
            
            b[ĩ] += q_r_n
            for s = 1:3
                j̃ = EToV[n, s]

                b_r, b_s = abc[r, 2], abc[s, 2]
                c_r, c_s = abc[r, 3], abc[s, 3]
                
                k_rs = 1 / (4 * abs(Δ)) * (lam1 * b_r * b_s + lam2 * c_r * c_s)
                A[ĩ, j̃] += k_rs
            end 
        end 
    end

    return A, b 
end 

println("CASE 1")
x0, y0 = 0, 0 
L1, L2 = 1, 1
noelms1, noelms2 = 4, 3 
qt = zeros(noelms1 * noelms2 * 2)
lam1 = 1 
lam2 = 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly_matrix(VX, VY, EToV, lam1, lam2, qt)
display(A)
display(diag(A))
#display(solve)

# function test_q(x,y)
#     return -6*x + 2*y - 2
# end

# qt_space = test_q.(VX,VY)