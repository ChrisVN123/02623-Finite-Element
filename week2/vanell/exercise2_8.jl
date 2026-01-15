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
include("exercise2_6.jl")

println("Exercise 2.8 a):")
x0, y0 = 0,0
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
lam1, lam2 = 1, 1 

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

u(x, y) = cos(pi*x)*sin(pi*y)
u_x(x, y) = ForwardDiff.derivative(u(x,y),x)
u_y(x, y) = ForwardDiff.derivative(u(x,y),x)
f(x, y) = u(x, y)
q̃(x, y) = 0 
q_l(x, y) = u_x(x, y)
q_b(x, y) = u_y(x, y)