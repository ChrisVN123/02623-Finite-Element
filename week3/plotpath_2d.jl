using CairoMakie

function plotting_path(EToV, solution, cols::Int, start_node::Int;
                       filename::Union{Nothing,String} = "maze_path.png",
                       tol::Float64 = 1e-12,
                       draw_walls::Bool = true,
                       wall_lw::Float64 = 6.0)

    u = collect(solution)
    n = length(u)

    # --- Build adjacency (undirected) ---
    adj = [Int[] for _ in 1:n]
    for e in EToV
        a, b = e[1], e[2]
        push!(adj[a], b)
        push!(adj[b], a)
    end

    # --- Greedy downhill path ---
    function greedy_downhill_path(adj, u; start::Int, tol::Float64)
        path = [start]
        current = start
        visited = Set([start])

        while true
            candidates = [v for v in adj[current] if (u[v] < u[current] - tol) && !(v in visited)]
            isempty(candidates) && break
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

    # --- Edge set for fast "is corridor open?" checks ---
    edgeset = Set{Tuple{Int,Int}}()
    for e in EToV
        a, b = e[1], e[2]
        a > b && ((a, b) = (b, a))
        push!(edgeset, (a, b))
    end
    has_edge(i, j) = (min(i, j), max(i, j)) in edgeset

    # --- Plot ---
    fig = Figure(size=(1000, 700))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xticks = 1:cols,
        yticks = 1:rows,
        title = "Greedy downhill path (min-u neighbor each step)"
    )

    # light grid lines (optional; keep if you like)
    for x in 1:cols
        vlines!(ax, x, color=:gray90)
    end
    for y in 1:rows
        hlines!(ax, y, color=:gray90)
    end

    # ---- Walls (draw missing neighbor edges as thick black walls) ----
    if draw_walls
        wall_color = :blue

        # Outer border
        lines!(ax, [0.5, cols + 0.5], [rows + 0.5, rows + 0.5], color=wall_color, linewidth=wall_lw) # top
        lines!(ax, [0.5, cols + 0.5], [0.5, 0.5],           color=wall_color, linewidth=wall_lw) # bottom
        lines!(ax, [0.5, 0.5],           [0.5, rows + 0.5], color=wall_color, linewidth=wall_lw) # left
        lines!(ax, [cols + 0.5, cols + 0.5], [0.5, rows + 0.5], color=wall_color, linewidth=wall_lw) # right

        # Interior walls between horizontal neighbors (i, i+1)
        for i in 1:n
            # right neighbor exists if not at end of row
            if xs[i] < cols
                j = i + 1
                if j <= n && !has_edge(i, j)
                    xwall = xs[i] + 0.5
                    y0 = ys[i] - 0.5
                    y1 = ys[i] + 0.5
                    lines!(ax, [xwall, xwall], [y0, y1], color=wall_color, linewidth=wall_lw)
                end
            end
        end

        # Interior walls between vertical neighbors (i, i+cols)
        for i in 1:(n - cols)
            j = i + cols
            if j <= n && !has_edge(i, j)
                ywall = (ys[i] + ys[j]) / 2   # should be ys[i] - 0.5 with your coordinate convention
                x0 = xs[i] - 0.5
                x1 = xs[i] + 0.5
                lines!(ax, [x0, x1], [ywall, ywall], color=wall_color, linewidth=wall_lw)
            end
        end
    end

    # all edges (light) - corridors (optional; comment out if you only want walls)
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

    # (REMOVED) labels: node id and u-value
    # text!(...)  <-- gone

    CairoMakie.xlims!(ax, 0.5, cols + 0.5)
    CairoMakie.ylims!(ax, 0.5, rows + 0.5)

    if filename !== nothing
        save(filename, fig)
    end

    return fig, path
end