using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")
include("exercise2_3.jl")

function calculate_bnodes(noelms1, noelms2)
    """
    Should work, please confirm 
    """
    bnodes = []
    for i ∈ 0:noelms1, j ∈ 0:noelms2
        if (i == 0 || i == noelms1 || j == 0 || j == noelms2)
            # push!(bnodes, (i, j, j+1 + i * (noelms1 + 1)))
            push!(bnodes, j+1 + i * (noelms2 + 1))
        end 
    end 
    return bnodes
end

function dirbc(bnodes, f, A, b)
    """
    Algorithm 6 
    It seems that Γ₂ = Γ 

    bnodes = global node numbers of the nodes on Γ₂ 
    """

    M = size(A)[1]
    for i ∈ bnodes 
        A[i, i] = 1 
        b[i] = f[i] 
        for j ∈ 1:M
            if j != i 
                A[i, j] = 0 
                if !(j ∈ bnodes)
                    b[j] += -A[j, i] * f[i]
                    A[j, i] = 0
                end 
            end 
        end 
    end 

    return A, b 
end 