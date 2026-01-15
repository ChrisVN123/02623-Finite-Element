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
    Make sure that û[j] matches u(xj, yj)! 
    """
    return maximum(
        abs.(û - u.(VX, VY))
    )
end


println("EXERCISE 2.7")
println("CASE 1:")
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
lam1, lam2 = 1, 1 

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

u(x, y) = 3 * x + 5 * y - 7
u_x(x, y) = 3 
u_y(x, y) = 5
f(x, y) = u(x, y)
q̃(x, y) = 0 
q_l(x, y) = u_x(x, y)
q_b(x, y) = u_y(x, y)

# GLOBAL ASSEMBLY 
A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
bnodes = calculate_bnodes(noelms1, noelms2)
beds = ConstructBeds(noelms1, noelms2, EToV)

Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

# IMPOSING NEUMANN 
b_l = neubc(VX, VY, EToV, Γ1[1:noelms2,:], q_l, b)
b_b = neubc(VX, VY, EToV, Γ1[noelms1:end,:], q_b, b_l) 

# IMPOSING DIRICHLET 
bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)

# PLOTS
û = A \ b
println("û_j in 2D:")
display(reshape(û,(noelms1+1, noelms2+1)))
fig = plot_fem_solution3d(VX,VY,EToV,û)
save(joinpath(pwd(), "plots/2.7_case1_results.png"), fig)

err_2_7_1 = compute_error(VX, VY, û, u)
println("E = $err_2_7_1")

VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
EToV_analytical = conelmtab(1_000, 1_000)
fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.7_case1_analytical.png"), fig_analytic)

println("CASE 2:")
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 32, 32
lam1, lam2 = 1, 1 

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)

EToV = conelmtab(noelms1, noelms2)

u(x, y) = sin(x)*sin(y)
u_x(x, y) = cos(x)*sin(y)
u_y(x, y) = sin(x)*cos(y)
f(x, y) = u(x, y)
q̃(x, y) = 2 * sin(x) * sin(y)
q_l(x, y) = u_x(x, y)
q_b(x, y) = u_y(x, y)

# GLOBAL ASSEMBLY 
A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
bnodes = calculate_bnodes(noelms1, noelms2)
beds = ConstructBeds(noelms1, noelms2, EToV)
Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

# IMPOSING NEUMANN 
b_l = neubc(VX, VY, EToV, Γ1[1:noelms2,:], q_l, b)
b_b = neubc(VX, VY, EToV, Γ1[noelms1:end,:], q_b, b_l) 

# IMPOSING DIRICHLET 
bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)

û = A \ b

err_2_7_2 = compute_error(VX,VY, û, u)
println("E = $err_2_7_2")

fig = plot_fem_solution3d(VX,VY,EToV,û)
save(joinpath(pwd(), "plots/2.7_case2_results.png"), fig)

VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
EToV_analytical = conelmtab(1_000, 1_000)
fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.7_case2_analytical.png"), fig_analytic)