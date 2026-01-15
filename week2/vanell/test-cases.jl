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
include("exercise2_5.jl")

println("EXERCISE 2.1")
println("a) CASE 1:")
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
println("VX = $VX")
println("VY = $VY")

println("b) CASE 1:")
display(conelmtab(4, 3))

println("\nEXERCISE 2.2")
println("CASE 2.2a:")
n = 4 # or 9? 
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
d, abc = basfun(n, VX, VY, EToV)

println("Δ = $d")
println("abc = ")
display(abc)

x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

println("CASE 2.2b:")
for k = 1:3 
    println("k = $k")
    n_vec = outernormal(n, k, VX, VY, EToV) # n = 4 
    display(n_vec)
end 

println("\nEXERCISE 2.3")
println("CASE 1:")
x0, y0 = 0, 0 
L1, L2 = 1, 1
noelms1, noelms2 = 4, 3 
qt = zeros(noelms1 * noelms2 * 2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)
display(A)

function spdiags(A)
    """
    A is square 
    Should work!
    """
    M = size(A)[1]
    B = zeros(M, M)
    d = []
    
    count = 1
    for k ∈ [-M:0; 1:M]
        kth_diag = diag(A, k)
        if !iszero(kth_diag)
            l = length(kth_diag)
            if k <= 0 
                B[1:l, count] = kth_diag 
            else 
                B[M-l+1:end, count] = kth_diag
            end 
            push!(d, k)
            count += 1
        end 
    end 

    return B[:, 1:count-1], d 
end

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)


println("CASE 2:")
q(x, y) = -6 * x + 2 * y - 2
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
qt = q.(VX, VY)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)

println("\nEXERCISE 2.4")
println("CASE 1:")
x0, y0 = 0, 0 
L1, L2 = 1, 1
noelms1, noelms2 = 4, 3
qt = zeros(noelms1 * noelms2 * 2)
f = ones(noelms1 * noelms2 * 2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
display(EToV)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)
bnodes = calculate_bnodes(noelms1, noelms2)
A, b = dirbc(bnodes, f, A, b)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)

println("Array(A)[1:13, 1:13] = ")
display(Array(A)[1:13, 1:13])

println("CASE 2:")
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
q(x, y) = -6 * x + 2 * y - 2 
f24_2(x, y) = x^3 - x^2 * y + y^2 - 1
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly(VX, VY, EToV, lam1, lam2, q.(VX, VY))
bnodes = calculate_bnodes(noelms1, noelms2)
A, b = dirbc(bnodes, f24_2.(VX, VY), A, b)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)

println("\nEXERCISE 2.5")
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


