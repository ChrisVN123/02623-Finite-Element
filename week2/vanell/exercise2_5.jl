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

"""
We consider the BVP 
u_{x x} + u_{y y} = -q̃(x, y), 
u = f(x, y). 
"""

function compute_error(VX, VY, û, u)
    """ 
    u is a function! 
    """
    return maximum(
        abs.(û - u.(VX, VY))
    )
end
