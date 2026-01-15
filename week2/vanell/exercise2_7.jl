using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")
include("exercise2_3.jl")
include("exercise2_4.jl")
include("plot.jl")
include("exercise2_6.jl")

"""
Γ₁ = left and bottom edges 
Γ₂ = right and top edges 

u_xx + u_yy = -q̃(x, y), (x, y) ∈ Ω

un = -q(x, y), (x, y) ∈ Γ₁ 
u = f(x, y),   (x, y) ∈ Γ₂
"""

function boundary_edges(noelms1, noelms2, beds)
    Γ1 = Matrix{Int}(undef, noelms1+noelms2, 2)
    Γ2 = Matrix{Int}(undef, noelms1+noelms2, 2)

    k1 = 1
    k2 = 1
    for (i, v) in enumerate(beds[:, 1])
        if v % 2 == 0 
            Γ1[k1, 1] = v
            Γ1[k1, 2] = beds[i, 2]
            k1 += 1
        elseif v % 2 == 1 
            Γ2[k2, 1] = v
            Γ2[k2, 2] = beds[i, 2]
            k2 += 1
        end 
    end 
    
    return Γ1, Γ2
end

function dirichlet_bound(noelms1, noelms2)   
    bnodes_dirbc = []
    for i in 0:(noelms1)
        push!(bnodes_dirbc, 1 + i * (noelms2+1))
    end 
    bnodes_dirbc = [bnodes_dirbc; collect(last(bnodes_dirbc)+1:last(bnodes_dirbc)+noelms2)]
    
    return bnodes_dirbc
end

function compute_error(VX, VY, û, u)
    """ 
    u is a function! 
    Estimating the error between the solution and the approximation at the nodes
    """
    return maximum(
        abs.(û - u.(VX, VY))
    )
end
