using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

function element_mesh_matrix(x0, y0, L1, L2, noelms1, noelms2)
    dx = L1 / noelms1
    dy = L2 / noelms2

    coords = Matrix{NTuple{3, Float64}}(undef, noelms2+1, noelms1+1)
    id = 1

    for j in 0:noelms1               # x index
        for i in 0:noelms2           # y index (top to bottom)
            x = x0 + j*dx
            y = y0 + L2 - i*dy
            coords[i+1, j+1] = (x, y, id)
            id += 1
        end
    end

    return coords
end

ptsM = element_mesh_matrix(0, 0, 4, 3, 4, 3)
# @show size(ptsM)   # (3, 3)
display(ptsM)


function EToV(nx::Int, ny::Int)
    @assert nx ≥ 2 && ny ≥ 2 
    nx += 1
    ny += 1
    num_tris = 2 * (nx - 1) * (ny - 1)
    etoV = Matrix{Int}(undef, num_tris, 4)

    e = 1
    tl = 1 

    for col in 1:(nx - 1)            
        tl = 1 + (col - 1) * ny   
        for row in 1:(ny - 1)          
            tr = tl + ny
            bl = tl + 1
            br = tl + ny + 1
            
            #tri1
            etoV[e, 1] = e
            etoV[e, 2] = tr
            etoV[e, 3] = tl
            etoV[e, 4] = br
            e += 1

            #tri2
            etoV[e, 1] = e
            etoV[e, 2] = bl
            etoV[e, 3] = br
            etoV[e, 4] = tl
            e += 1

            tl += 1  #move one node down within the column
        end
    end

    return etoV
end

display(EToV(4,3))

function delta_abc(n, VX, VY, EToV, coordinates)
    # values[e, k] holds (x_k, y_k) for element e, k = 1..3
    values = Matrix{NTuple{2,Float64}}(undef, n, 3)

    # --- gather triangle vertex coordinates (using columns 2:4 of EToV) ---
    for e in 1:n, j in 2:4
        vid = EToV[e, j]  # vertex id
        idx = findfirst(t -> t[end] == vid, coordinates)
        idx === nothing && error("Vertex id $vid not found in coordinates.")
        row, col = Tuple(idx)
        x, y = coordinates[row, col][1:2]
        values[e, j-1] = (x, y)
    end

    # abc[e, i, :] = [ã_i, b̃_i, c̃_i], i = 1..3
    abc   = Array{Float64}(undef, n, 3, 3)
    delta = Vector{Float64}(undef, n)

    for e in 1:n
        x1, y1 = values[e, 1]
        x2, y2 = values[e, 2]
        x3, y3 = values[e, 3]

        # signed area (Delta)
        delta[e] = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        # (i,j,k) = (1,2,3)
        abc[e, 1, 1] = x2*y3 - x3*y2   # ã1
        abc[e, 1, 2] = y2 - y3         # b̃1
        abc[e, 1, 3] = x3 - x2         # c̃1

        # (i,j,k) = (2,3,1)
        abc[e, 2, 1] = x3*y1 - x1*y3   # ã2
        abc[e, 2, 2] = y3 - y1         # b̃2
        abc[e, 2, 3] = x1 - x3         # c̃2

        # (i,j,k) = (3,1,2)
        abc[e, 3, 1] = x1*y2 - x2*y1   # ã3
        abc[e, 3, 2] = y1 - y2         # b̃3
        abc[e, 3, 3] = x2 - x1         # c̃3
    end

    display(delta)
    display(abc)
    return abc, delta
end


etoV = EToV(4, 3)
n = size(etoV, 1)
abc, delta = delta_abc(n, nothing, nothing, etoV, ptsM)


delta_abc(24, 1,1, EToV(4,3),ptsM)






function plot_mesh(ptsM, etoV; show_ids=true)
    ny, nx = size(ptsM)

    # 1) Build id -> (x,y)
    id2xy = Dict{Int, Tuple{Float64,Float64}}()
    xs = Float64[]
    ys = Float64[]
    ids = Int[]

    for I in CartesianIndices(ptsM)
        x, y, idf = ptsM[I]
        id = Int(idf)
        id2xy[id] = (x, y)
        push!(xs, x); push!(ys, y); push!(ids, id)
    end

    # 2) Collect unique edges from triangles
    # Each triangle has edges (a,b), (b,c), (c,a)
    edges = Set{Tuple{Int,Int}}()
    @inbounds for r in 1:size(etoV, 1)
        a = etoV[r, 2]; b = etoV[r, 3]; c = etoV[r, 4]
        push!(edges, (min(a,b), max(a,b)))
        push!(edges, (min(b,c), max(b,c)))
        push!(edges, (min(c,a), max(c,a)))
    end

    # 3) Plot nodes
    p = scatter(xs, ys;
        legend=false,
        aspect_ratio=:equal,
        grid=true,
        xlabel="x", ylabel="y",
        title="FEM mesh"
    )

    # 4) Draw edges
    for (u, v) in edges
        (x1, y1) = id2xy[u]
        (x2, y2) = id2xy[v]
        plot!(p, [x1, x2], [y1, y2]; label=false)
    end

    # 5) Optional node labels
    if show_ids
        for (x, y, id) in zip(xs, ys, ids)
            annotate!(p, x, y, text(string(id), 9, :left))
        end
    end

    return savefig("ex2_2.png")
end
etoV = EToV(4, 3)          # nx=5, ny=4 (nodes, not elements)
plot_mesh(ptsM, etoV)

