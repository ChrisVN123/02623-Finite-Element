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
f(x, y) = 3 * x + 5 * y - 7
u_x(x, y) = 3
u_y(x, y) = 5
q(x, y) = 0

A, b = assembly(VX, VY, EToV, lam1, lam2, q.(VX, VY))
println("size A, b = $(size(A)), $(size(b))")
bnodes = calculate_bnodes(noelms1, noelms2)

CB = ConstructBeds(VX, VY, EToV, "")
println("CB = ")
display(CB)

Γ1 = Matrix{Int}(undef, noelms1+noelms2, 2)
Γ2 = Matrix{Int}(undef, noelms1+noelms2, 2)

k1 = 1
k2 = 1

for (i, v) in enumerate(CB[:,1])
    println("i, v = $i, $v")
    if v % 2 == 0 
        println(v, v % 2) 
        # push!(Γ1, v)
        Γ1[k1, 1] = v
        Γ1[k1, 2] = CB[i, 2]
        global k1 += 1
    elseif v % 2 == 1 
        # push!(Γ2, v)
        Γ2[k2, 1] = v
        Γ2[k2, 2] = CB[i, 2]
        global k2 += 1
    end 
end 

println(Γ1)
println(Γ2)

q_l(x, y) = -u_x(x, y)
q_b(x, y) = -u_y(x, y)

b_l = neubc(VX, VY, EToV, Γ1, q_l, b)
b_b = neubc(VX, VY, EToV, Γ1, q_b, b_l) # maybe! 

println("size A, b = $(size(A)), $(size(b_b))")
println("Γ2 = $Γ2")

bnodes_dirbc = []
for i in 0:(noelms1)
    push!(bnodes_dirbc, 1 + i * (noelms2+1))
end 
# println(collect(last(bnodes_dirbc)+1:last(bnodes_dirbc)+noelms2))
bnodes_dirbc = [bnodes_dirbc; collect(last(bnodes_dirbc)+1:last(bnodes_dirbc)+noelms2)]
println("TEST = $bnodes_dirbc")

A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)
û = A \ b
fig = plot_fem_solution3d(VX,VY,EToV,û)
save(joinpath(pwd(), "plots/2.7_results.png"), fig)

println(û)

@time VX_analytical, VY_analytical = xy(x0, y0, L1, L2, 1_000, 1_000)
@time EToV_analytical = conelmtab(1_000, 1_000)
@time fig_analytic = plot_fem_solution3d(VX_analytical, VY_analytical, EToV_analytical, u.(VX_analytical, VY_analytical))
save(joinpath(pwd(), "plots/2.7_case1_analytical.png"), fig_analytic)