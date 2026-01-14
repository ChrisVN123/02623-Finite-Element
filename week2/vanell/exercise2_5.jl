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
    Estimating the error between the solution and the approximation at the nodes
    Make sure that û[j] matches u(xj, yj)! 
    """
    return maximum(
        abs.(û - u.(VX, VY))
    )
end

println("CASE 1:")
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
lam1, lam2 = 1, 1 

u(x, y) = x^3 - x^2 * y + y^2 - 1 
q̃(x, y) = -6 * x + 2 * y - 2 
f(x, y) = u(x, y)

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
bnodes = calculate_bnodes(noelms1, noelms2)
A, b = dirbc(bnodes, f.(VX, VY), A, b)
û = A \ b 

println("û = $û")
println("E = $(compute_error(VX, VY, û, u))")

fig_solution = plot_fem_solution3d(VX, VY, EToV, û)
save(joinpath(pwd(), "plots/2.5_case1_solution.png"), fig_solution)

VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
EToV_analytical = conelmtab(1_000, 1_000)
fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.5_case1_analytical.png"), fig_analytic)

println("CASE 2:")
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 

u(x, y) = x^2 * y^2 
q̃(x, y) = -2 * x^2 - 2 * y^2
f(x, y) = u(x, y)

E_arr = []
h_arr = []
û_arr = []
for p in 1:6
    println("p = $p")
    local noelms1, noelms2 = 2^p, 2^p 
    
    push!(h_arr, √((L2/noelms2)^2+(L1/noelms1)^2))
    local VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    local EToV = conelmtab(noelms1, noelms2)
    local A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
    local bnodes = calculate_bnodes(noelms1, noelms2)
    local A, b = dirbc(bnodes, f.(VX, VY), A, b)
    local û = A \ b 

    push!(û_arr, (p, û))
    
    fig = plot_fem_solution3d(VX, VY, EToV, u.(VX, VY))
    save(joinpath(pwd(), "plots/2.5_case2_p=$(p).png"), fig)
    push!(E_arr, compute_error(VX, VY, û, u))
end 

# println("û_arr = $û_arr") # prints a LOT 
println("E_arr = $E_arr")

Plots.plot(log.(h_arr), log.(E_arr), title="E(p)", label="Error")
hs = LinRange(log(minimum(h_arr)), log(maximum(h_arr)), 1000)
ys = 0.997591 .+ 2 * hs
Plots.plot!(hs, ys, label = "Theoretical result")
savefig("plots/2.5_case2_log_error.png")

println("polynomial fit: \n$(fit(log.(h_arr), log.(E_arr), 1))")

VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
EToV_analytical = conelmtab(1_000, 1_000)
fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.5_case2_analytical.png"), fig_analytic)


