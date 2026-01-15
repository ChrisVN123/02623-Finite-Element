using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials
using Pkg

include("exercise2_1.jl")
include("exercise2_2.jl")
include("exercise2_3.jl")
include("exercise2_4.jl")
include("plot.jl")
#include("exercise2_5.jl")
include("exercise2_6.jl")

"""
Γ₁ = left and bottom edges 
Γ₂ = right and top edges 

u_xx + u_yy = -q̃(x, y), (x, y) ∈ Ω

un = -q(x, y), (x, y) ∈ Γ₁ 
u = f(x, y),   (x, y) ∈ Γ₂
"""

println("EXERCISE 2.7")
println("CASE 1:")
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
lam1, lam2 = 1, 1 

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

u(x, y) = 3 * x + 5 * y - 7
f(x, y) = u(x, y)
q̃(x, y) = 0 
u_x(x, y) = 3
u_y(x, y) = 5

A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))

function boundary_edge_pairs(VX, VY, EToV, x0, y0, L1, L2; tol=1e-10)
    left = Tuple{Int, Int}[]
    bottom = Tuple{Int, Int}[]
    right = Tuple{Int, Int}[]
    top = Tuple{Int, Int}[]
    x_left = x0
    x_right = x0 + L1
    y_bottom = y0
    y_top = y0 + L2

    for n in 1:size(EToV, 1)
        for r in 1:3
            s = r == 1 ? 2 : r == 2 ? 3 : 1
            i, j = EToV[n, r], EToV[n, s]
            xi, yi = VX[i], VY[i]
            xj, yj = VX[j], VY[j]

            if abs(xi - x_left) <= tol && abs(xj - x_left) <= tol
                push!(left, (n, r))
            elseif abs(xi - x_right) <= tol && abs(xj - x_right) <= tol
                push!(right, (n, r))
            elseif abs(yi - y_bottom) <= tol && abs(yj - y_bottom) <= tol
                push!(bottom, (n, r))
            elseif abs(yi - y_top) <= tol && abs(yj - y_top) <= tol
                push!(top, (n, r))
            end
        end
    end

    return left, bottom, right, top
end

function pairs_to_matrix(pairs)
    out = Matrix{Int}(undef, length(pairs), 2)
    for (idx, pair) in enumerate(pairs)
        out[idx, 1] = pair[1]
        out[idx, 2] = pair[2]
    end
    return out
end

left_pairs, bottom_pairs, _, _ = boundary_edge_pairs(VX, VY, EToV, x0, y0, L1, L2)
beds_left = pairs_to_matrix(left_pairs)
beds_bottom = pairs_to_matrix(bottom_pairs)

q_l(x, y) = u_x(x, y)
q_b(x, y) = u_y(x, y)

println("b = $b")
b_l = neubc(VX, VY, EToV, beds_left, q_l, b)
println("b_l = $b_l")
b_b = neubc(VX, VY, EToV, beds_bottom, q_b, b_l)

println("b_b = $b_b")

bnodes_dirbc = []
for i in 0:(noelms1)
    push!(bnodes_dirbc, 1 + i * (noelms2+1))
end 
bnodes_dirbc = [bnodes_dirbc; collect(last(bnodes_dirbc)+1:last(bnodes_dirbc)+noelms2)]
println("b_dbc = $bnodes_dirbc")
A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)
û = A \ b
fig = plot_fem_solution3d(VX,VY,EToV,û)
save(joinpath(pwd(), "plots/2.7_results_new.png"), fig)
print(b)

@time VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
@time EToV_analytical = conelmtab(1_000, 1_000)
@time fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.7_case1_analytical.png"), fig_analytic)
