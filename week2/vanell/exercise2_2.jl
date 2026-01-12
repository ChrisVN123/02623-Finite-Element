using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials
include("exercise2_1.jl")

function basfun(n, VX, VY, EToV)
    abc = Matrix{Float64}(undef, 3, 3)
    
    i, j, k = EToV[n, :]
    x1, y1 = VX[i], VY[i]
    x2, y2 = VX[j], VY[j]
    x3, y3 = VX[k], VY[k]

    Δ = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

    # (i,j,k) = (1, 2, 3)
    abc[1, 1] = x2*y3 - x3*y2   # ã₁
    abc[1, 2] = y2 - y3         # b̃₁
    abc[1, 3] = x3 - x2         # c̃₁

    # (i,j,k) = (2, 3, 1)
    abc[2, 1] = x3*y1 - x1*y3   # ã₂
    abc[2, 2] = y3 - y1         # b̃₂
    abc[2, 3] = x1 - x3         # c̃₂

    # (i,j,k) = (3, 1, 2)
    abc[3, 1] = x1*y2 - x2*y1   # ã₃ 
    abc[3, 2] = y1 - y2         # b̃₃
    abc[3, 3] = x2 - x1         # c̃₃
    return Δ, abc
end


# display(EToV)
println("CASE 2.2a:")
n = 4
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
d, abc = basfun(n, VX, VY, EToV)
display(d)
display(abc)

function outernormal(n, k, VX, VY, EToV)
    i, j, h = EToV[n, :]
    x1, x2 = nothing, nothing
    y1, y2 = nothing, nothing
    
    if k == 1 
        x1, y1 = VX[i], VY[i]
        x2, y2 = VX[j], VY[j]
    elseif k == 2 
        x1, y1 = VX[j], VY[j]
        x2, y2 = VX[h], VY[h]
    elseif k == 3 
        x1, y1 = VX[h], VY[h] 
        x2, y2 = VX[i], VY[i]
    end 
    
    t1 = x2 - x1
    t2 = y2 - y1
    normal = 1 / (t1^2 + t2^2)^0.5
    
    return [t2, -t1] * normal
end

x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

println("CASE 2.2b:")
for k = 1:3 
    println("k = $k")
    n_vec = outernormal(9, k, VX, VY, EToV)
    display(n_vec)
end 

# print("Running Element mesh creation")
# @time ptsM = xy(0.0, 0.0, 4.0, 3.0, 4, 3)
# print("Running EToV")
# @time etoV = EToV(4, 3)

# n = size(etoV, 1)
# abc, delta = delta_abc(n, nothing, nothing, etoV, ptsM)


#delta_abc(24, 1,1, EToV(4,3),ptsM)











# function plot_mesh(ptsM, etoV; show_ids=true)
#     ny, nx = size(ptsM)

#     # 1) Build id -> (x,y)
#     id2xy = Dict{Int, Tuple{Float64,Float64}}()
#     xs = Float64[]
#     ys = Float64[]
#     ids = Int[]

#     for I in CartesianIndices(ptsM)
#         x, y, idf = ptsM[I]
#         id = Int(idf)
#         id2xy[id] = (x, y)
#         push!(xs, x); push!(ys, y); push!(ids, id)
#     end

#     # 2) Collect unique edges from triangles
#     # Each triangle has edges (a,b), (b,c), (c,a)
#     edges = Set{Tuple{Int,Int}}()
#     @inbounds for r in 1:size(etoV, 1)
#         a = etoV[r, 2]; b = etoV[r, 3]; c = etoV[r, 4]
#         push!(edges, (min(a,b), max(a,b)))
#         push!(edges, (min(b,c), max(b,c)))
#         push!(edges, (min(c,a), max(c,a)))
#     end

#     # 3) Plot nodes
#     p = scatter(xs, ys;
#         legend=false,
#         aspect_ratio=:equal,
#         grid=true,
#         xlabel="x", ylabel="y",
#         title="FEM mesh"
#     )

#     # 4) Draw edges
#     for (u, v) in edges
#         (x1, y1) = id2xy[u]
#         (x2, y2) = id2xy[v]
#         plot!(p, [x1, x2], [y1, y2]; label=false)
#     end

#     # 5) Optional node labels
#     if show_ids
#         for (x, y, id) in zip(xs, ys, ids)
#             annotate!(p, x, y, text(string(id), 9, :left))
#         end
#     end

#     return savefig("ex2_2.png")
# end
# etoV = EToV(4, 3)          # nx=5, ny=4 (nodes, not elements)
# plot_mesh(ptsM, etoV)

