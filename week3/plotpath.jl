using CairoMakie

function plotting_path(EToV, solution, cols::Int, start_node::Int;
                       filename::Union{Nothing,String} = "images/maze_path_old.png",
                       tol::Float64 = 1e-12)

    u = collect(solution)               # ensure indexable Vector
    n = length(u)

    # --- Build adjacency (undirected) ---
    adj = [Int[] for _ in 1:n]
    for e in EToV
        a, b = e[1], e[2]               # works for [a,b] vectors
        push!(adj[a], b)
        push!(adj[b], a)
    end

    # --- Greedy downhill path ---
    function greedy_downhill_path(adj, u; start::Int, tol::Float64)
        path = [start]
        current = start
        visited = Set([start])  # prevents cycling

        while true
            candidates = [v for v in adj[current] if (u[v] < u[current] - tol) && !(v in visited)]
            isempty(candidates) && break

            # pick neighbor with minimum u
            nextv = candidates[argmin(u[candidates])]
            push!(path, nextv)
            push!(visited, nextv)
            current = nextv
        end
        return path
    end

    path = greedy_downhill_path(adj, u; start=start_node, tol=tol)

    # --- Grid coordinates (row-major) ---
    rows = cld(n, cols)
    xs = [((i - 1) % cols) + 1 for i in 1:n]
    ys = [rows - ((i - 1) ÷ cols) for i in 1:n]  # row 1 at top

    # --- Plot ---
    fig = Figure(size=(1000, 700))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xticks = 1:cols,
        yticks = 1:rows,
        title = "Greedy downhill path (min-u neighbor each step)"
    )

    # grid lines
    for x in 1:cols
        vlines!(ax, x, color=:gray90)
    end
    for y in 1:rows
        hlines!(ax, y, color=:gray90)
    end

    # all edges (light)
    for e in EToV
        a, b = e[1], e[2]
        lines!(ax, [xs[a], xs[b]], [ys[a], ys[b]], linewidth=2, color=:gray70)
    end

    # path edges (highlight)
    if length(path) > 1
        for k in 1:(length(path)-1)
            a, b = path[k], path[k+1]
            lines!(ax, [xs[a], xs[b]], [ys[a], ys[b]], linewidth=6, color=:red)
        end
    end

    # nodes
    CairoMakie.scatter!(ax, xs, ys, markersize=14, color=:dodgerblue)

    # labels: node id and u-value
    labels = ["$i\n$(round(u[i]; digits=3))" for i in 1:n]
    text!(ax, labels, position=Point2f.(xs .+ 0.12, ys .+ 0.12), fontsize=12)

    CairoMakie.xlims!(ax, 0.5, cols + 0.5)
    CairoMakie.ylims!(ax, 0.5, rows + 0.5)

    # save (optional)
    if filename !== nothing
        save(filename, fig)
    end

    return fig, path
end


using CairoMakie

function plot_bfs_edges(EToV, bfs0, cols::Int; filename::String="maze_bfs.png", label_steps::Bool=true)
    # bfs0 entries are 0-based indices into EToV rows -> convert to 1-based
    bfs_edges = bfs0 .+ 1

    # Number of nodes inferred from EToV
    n = maximum(vcat(EToV...))
    rows = cld(n, cols)

    # Grid coordinates (row-major)
    xs = [((i - 1) % cols) + 1 for i in 1:n]
    ys = [rows - ((i - 1) ÷ cols) for i in 1:n]  # row 1 at top

    fig = Figure(size=(1200, 800))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xticks = 1:cols,
        yticks = 1:rows,
        title = "BFS edge order overlay (bfs0 indexes rows of EToV)"
    )

    # Grid
    for x in 1:cols
        vlines!(ax, x, color=:gray90)
    end
    for y in 1:rows
        hlines!(ax, y, color=:gray90)
    end

    # All edges (light)
    for e in EToV
        a, b = e[1], e[2]
        lines!(ax, [xs[a], xs[b]], [ys[a], ys[b]], linewidth=2, color=:gray70)
    end

    # Nodes (light)
    CairoMakie.scatter!(ax, xs, ys, markersize=8, color=:gray40)

    # Overlay BFS-used edges (in order)
    for (k, ei) in enumerate(bfs_edges)
        a, b = EToV[ei][1], EToV[ei][2]
        lines!(ax, [xs[a], xs[b]], [ys[a], ys[b]], linewidth=6, color=:red)

        if label_steps
            mx = (xs[a] + xs[b]) / 2
            my = (ys[a] + ys[b]) / 2
            text!(ax, string(k-1), position=Point2f(mx + 0.08, my + 0.08), fontsize=10)
        end
    end

    # Mark approximate start/end (based on first/last BFS edge)
    if !isempty(bfs_edges)
        a0, b0 = EToV[first(bfs_edges)][1], EToV[first(bfs_edges)][2]
        aL, bL = EToV[last(bfs_edges)][1],  EToV[last(bfs_edges)][2]

        CairoMakie.scatter!(ax, [xs[a0]], [ys[a0]], markersize=18, color=:green)   # start-ish
        CairoMakie.scatter!(ax, [xs[bL]], [ys[bL]], markersize=18, color=:purple)  # end-ish
    end

    CairoMakie.xlims!(ax, 0.5, cols + 0.5)
    CairoMakie.ylims!(ax, 0.5, rows + 0.5)

    save(filename, fig)
    return fig
end
