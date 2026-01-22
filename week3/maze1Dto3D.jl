using Pkg
using Random

using SparseArrays, LinearAlgebra
using GLMakie

using SparseArrays

function assemble_K_graph(EToV, start, stop, cval, dval)
    # Robust node count (works for tuples or 2-vectors)
    num_nodes = maximum(Iterators.flatten(EToV))

    I = Int[]
    J = Int[]
    V = Float64[]
    b = zeros(num_nodes)

    for (i, j) in EToV
        push!(I, i); push!(J, i); push!(V,  1.0)
        push!(I, j); push!(J, j); push!(V,  1.0)
        push!(I, i); push!(J, j); push!(V, -1.0)
        push!(I, j); push!(J, i); push!(V, -1.0)
    end

    A = sparse(I, J, V, num_nodes, num_nodes)

    # Dirichlet RHS
    b[start] = cval
    b[stop]  = dval

    # Dirichlet BC
    A[start, :] .= 0
    A[stop,  :] .= 0
    A[start, start] = 1.0
    A[stop,  stop]  = 1.0

    dropzeros!(A)
    return A, b
end

EToV  = [[1,2],[2,3],[2,4],[4,5]]
nodes = [
    0.0  0.0  0.0;   # 1
    1.0  0.0  0.0;   # 2
    2.0  0.0  0.0;   # 3
    1.0 -1.0  0.0;   # 4
    1.0 -1.0 -1.0    # 5 (unused by this EToV)
]

# start, stop = 1, 5
# c, d  = 1, 0
# A, b = assemble_K_graph(EToV, start, stop, c, d)
# display(A)
# display(b)
# u = A \ b

# display(u)

function adjacency(EToV, n_nodes)
    adj = [Int[] for _ in 1:n_nodes]
    for (i,j) in EToV
        push!(adj[i], j)
        push!(adj[j], i)
    end
    return adj
end


function descent_path_robust(EToV, u; start, stop)
    n_nodes = length(u)
    adj = adjacency(EToV, n_nodes)

    path = [start]
    current = start
    visited = Set([start])

    for _ in 1:(10n_nodes)
        current == stop && break
        nbrs = adj[current]
        isempty(nbrs) && break

        # Prefer unvisited neighbors; fall back to all neighbors
        candset = [v for v in nbrs if !(v in visited)]
        if isempty(candset)
            candset = nbrs
        end

        # Pick neighbor with minimum u (ties: first)
        vals = u[candset]
        next = candset[argmin(vals)]

        push!(path, next)
        current = next
        push!(visited, current)
    end

    return path
end

# ---------------------------
# Plotting
# ---------------------------
function plot_network_3d(EToV, nodes, u; start=1, stop=2, filename="3D_test", path=Int[])
    fig = Figure(size=(900, 700))
    ax  = Axis3(fig[1, 1], xlabel="x", ylabel="y", zlabel="z", aspect=:data)

    umin, umax = minimum(u), maximum(u)

    # Draw edges, colored by u along the segment
    for (i, j) in EToV
        p1 = Point3f(nodes[i,1], nodes[i,2], nodes[i,3])
        p2 = Point3f(nodes[j,1], nodes[j,2], nodes[j,3])

        # Color gradient along line via endpoint values
        lines!(ax, [p1, p2];
            color = [u[i], u[j]],
            colormap = :viridis,
            colorrange = (umin, umax),
            linewidth = 4
        )
    end

    # Draw nodes colored by u
    pts = [Point3f(nodes[i,1], nodes[i,2], nodes[i,3]) for i in 1:size(nodes,1)]
    scatter!(ax, pts; markersize=16, color=u, colormap=:viridis, colorrange=(umin, umax))

    # Label nodes (optional)
    for i in 1:size(nodes,1)
        text!(ax, string(i), position=pts[i], align=(:left, :bottom), fontsize=14)
    end

    # Highlight start/stop
    scatter!(ax, [pts[start]]; markersize=26, color=:white)
    scatter!(ax, [pts[stop]];  markersize=26, color=:black)

    # Highlight path if provided
    if !isempty(path)
        ppts = [Point3f(nodes[i,1], nodes[i,2], nodes[i,3]) for i in path]
        lines!(ax, ppts; color=:red, linewidth=8)
    end

    Colorbar(fig[1, 2], colormap=:viridis, limits=(umin, umax), label="Potential u")


    filename = "$(filename)_3D_4views_cols.png"
    if filename !== nothing
        GLMakie.save(filename, fig)
    end
    return fig, path
end



# # ---------------------------
# 3D grid maze generator (perfect maze + optional loops)
# ---------------------------
# Map (x,y,z) -> node id (1-based)
idx(x, y, z, nx, ny) = x + (y-1)*nx + (z-1)*nx*ny

function generate_3d_maze(nx, ny, nz; seed=1234, loop_prob=0.10, spacing=(1.0, 1.0, 1.0))
    rng = MersenneTwister(seed)
    n_nodes = nx * ny * nz

    # Node coordinates
    nodes = Array{Float64}(undef, n_nodes, 3)
    for z in 1:nz, y in 1:ny, x in 1:nx
        id = idx(x,y,z,nx,ny)
        nodes[id,1] = (x-1) * spacing[1]
        nodes[id,2] = (y-1) * spacing[2]
        nodes[id,3] = (z-1) * spacing[3]
    end

    # Helper: valid neighbors in 3D grid
    dirs = ((1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1))
    function neighbors_of(x,y,z)
        nbrs = Tuple{Int,Int,Int}[]
        for (dx,dy,dz) in dirs
            xx, yy, zz = x+dx, y+dy, z+dz
            if 1 <= xx <= nx && 1 <= yy <= ny && 1 <= zz <= nz
                push!(nbrs, (xx,yy,zz))
            end
        end
        return nbrs
    end

    # DFS (recursive backtracker) to create a spanning tree (perfect maze)
    visited = falses(n_nodes)
    stack = Tuple{Int,Int,Int}[]
    start_cell = (1,1,1)
    push!(stack, start_cell)
    visited[idx(start_cell..., nx, ny)] = true

    edges = Set{Tuple{Int,Int}}()  # store undirected edges with (min,max)

    while !isempty(stack)
        (x,y,z) = stack[end]
        cur = idx(x,y,z,nx,ny)

        unvis = Tuple{Int,Int,Int}[]
        for (xx,yy,zz) in neighbors_of(x,y,z)
            nid = idx(xx,yy,zz,nx,ny)
            if !visited[nid]
                push!(unvis, (xx,yy,zz))
            end
        end

        if isempty(unvis)
            pop!(stack)
        else
            nxt = unvis[rand(rng, 1:length(unvis))]
            nxtid = idx(nxt..., nx, ny)

            a, b = min(cur, nxtid), max(cur, nxtid)
            push!(edges, (a,b))

            visited[nxtid] = true
            push!(stack, nxt)
        end
    end

    # Optionally add extra edges (loops) among adjacent grid neighbors
    for z in 1:nz, y in 1:ny, x in 1:nx
        cur = idx(x,y,z,nx,ny)
        for (xx,yy,zz) in neighbors_of(x,y,z)
            nxt = idx(xx,yy,zz,nx,ny)
            a, b = min(cur, nxt), max(cur, nxt)
            if !( (a,b) in edges ) && rand(rng) < loop_prob
                push!(edges, (a,b))
            end
        end
    end

    # Convert to EToV vector of pairs
    EToV = [(a,b) for (a,b) in edges]

    # Define endpoints
    start = idx(1,1,1,nx,ny)
    stop  = idx(nx,ny,nz,nx,ny)

    return EToV, nodes, start, stop
end

# ---------------------------
# Run a "more complex" test
# ---------------------------
nx, ny, nz = 4, 4, 4                 # increase for more complexity
c, d  = 1, 0
EToV, nodes, start, stop = generate_3d_maze(nx, ny, nz; seed=2026, loop_prob=0.12, spacing=(1.0, 1.0, 1.0))


A, b = assemble_K_graph(EToV, start, stop, c, d)


u = A \ b

path = descent_path_robust(EToV, u; start=start, stop=stop)

fig = plot_network_3d(EToV, nodes, u; start=start, stop=stop, path=path)
display(fig)

println("nodes = ", size(nodes,1), ", edges = ", length(EToV))
println("start = ", start, ", stop = ", stop)
println("path length = ", length(path), ", reached stop = ", (last(path) == stop))