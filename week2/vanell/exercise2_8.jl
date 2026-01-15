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
include("exercise2_7.jl")

println("EXERCISE 2.8")
println("a)")

println("b)")
function program_b(n, name)
    x0, y0 = 0, 0
    L1, L2 = 1, 1
    noelms1, noelms2 = n, n
    lam1, lam2 = 1, 1 

    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV = conelmtab(noelms1, noelms2)

    f(x, y) = cos(pi * x) * cos(pi * y) # u(x, y) = f(x, y), (x, y) ∈ Γ2
    u_x(x, y) = 0
    u_y(x, y) = 0
    q̃(x, y) = 2 * pi^2 * cos(pi * x) * cos(pi * y)
    q_l(x, y) = u_x(x, y)
    q_b(x, y) = u_y(x, y)

    # GLOBAL ASSEMBLY 
    A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
    bnodes = calculate_bnodes(noelms1, noelms2)
    beds = ConstructBeds(noelms1, noelms2, EToV)
    Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

    # IMPOSING NEUMANN 
    b_l = neubc(VX, VY, EToV, Γ1[1:noelms2, :], q_l, b)
    b_b = neubc(VX, VY, EToV, Γ1[noelms1:end, :], q_b, b_l) 

    # IMPOSING DIRICHLET 
    bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
    A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)

    # SOLUTION
    û = A \ b

    fig = plot_fem_solution3d(VX, VY, EToV, û)
    save(joinpath(pwd(), "plots/$name.png"), fig)
    return û[1+noelms2]
end
program_b(3, "2.8_b_results")

println("c)")
function program_c(n, name)
    x0, y0 = -1, -1
    L1, L2 = 2, 2
    noelms1, noelms2 = n, n
    lam1, lam2 = 1, 1 

    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV = conelmtab(noelms1, noelms2)

    f(x, y) = cos(pi * x) * cos(pi * y) # u(x, y) = f(x, y), (x, y) ∈ Γ2
    q̃(x, y) = 2 * pi^2 * cos(pi * x) * cos(pi * y)

    # GLOBAL ASSEMBLY 
    A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
    bnodes = calculate_bnodes(noelms1, noelms2)
    beds = ConstructBeds(noelms1, noelms2, EToV)
    Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

    # IMPOSING DIRICHLET 
    bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
    A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b)

    # SOLUTION
    û = A \ b

    fig = plot_fem_solution3d(VX, VY, EToV, û)
    save(joinpath(pwd(), "plots/$name.png"), fig)    
    
    idx_0 = Int(ceil(length(û) / 2))
    return VX, VY, û, û[idx_0]
end
program_c(6, "2.8_c_results")

println("d)")
p_arr = collect(1:5)
err_arr_c = []
u_c(x, y) = cos(pi * x) * cos(pi * y)

for p ∈ p_arr
    nb = 2^p 
    nc = 2^(p+1)

    û0_b = program_b(nb, "2.8_d_b_$p")
    local VX, VY, û_c, û0_c = program_c(nc, "2.8_d_c_$p")

    err = compute_error(VX, VY, û_c, u_c)
    push!(err_arr_c, err)

    println("p = $p")
    println("\tû0_b = $û0_b")
    println("\tû0_c = $û0_c")
end 

println("err_arr_c = $err_arr_c")

println("fit: $(fit(
    log.(2 .^ (p_arr .+ 1)), 
    log.(err_arr_c), 1)
)")