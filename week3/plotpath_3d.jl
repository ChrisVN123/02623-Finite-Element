using GLMakie
import Makie

function plotting_path_3d(EToV, solution, cols::Int, start_node::Int;
                          filename::Union{Nothing,String} = "maze",
                          tol::Float64 = 1e-12,
                          zscale::Float64 = 1.0)

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

    # --- Grid coordinates (row-major), but NOT flipped in y for 3D ---
    rows = cld(n, cols)
    xs = [((i - 1) % cols) + 1 for i in 1:n]
    ys = [((i - 1) ÷ cols) + 1 for i in 1:n]   # row 1 at "top" but now y=1 (front)

    # --- Surface Z(row, col) = u ---
    Z = fill(NaN, rows, cols)
    for i in 1:n
        Z[ys[i], xs[i]] = zscale * u[i]
    end
    xgrid = collect(1:cols)
    ygrid = collect(1:rows)

    fig = GLMakie.Figure(size=(1100, 800))
    ax = GLMakie.Axis3(fig[1, 1];
        title="3D maze potential (height = u)",
        xlabel="x", ylabel="y", zlabel="z = u",
        aspect = :data
    )

    Makie.surface!(ax, xgrid, ygrid, permutedims(Z))
    #Makie.scatter!(ax, xs, ys, zscale .* u; markersize=10)

    if length(path) > 1
        #path_pts = GLMakie.Point3f.([xs[i] for i in path],
        #                            [ys[i] for i in path],
        #                            zscale .* [u[i] for i in path])
        #Makie.lines!(ax, path_pts; linewidth=4, color=:red)
    end

    Makie.xlims!(ax, 0.5, cols + 0.5)
    Makie.ylims!(ax, 0.5, rows + 0.5)

    # Optional: choose a nicer default view
    # ax.azimuth[] = -pi/2
    # ax.elevation[] =  pi/6
    filename = "$(filename)_3D_cols_$(cols).png"
    if filename !== nothing
        Makie.save(filename, fig)
    end

    return fig, path
end

using GLMakie
import Makie

function plotting_path_3d_4views(EToV, solution, cols::Int, start_node::Int;
                                 filename::Union{Nothing,String} = "maze",
                                 tol::Float64 = 1e-12,
                                 zscale::Float64 = 1.0,
                                 markersize::Float64 = 10.0)

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

    # --- Grid coordinates (row-major), NOT flipped for 3D ---
    rows = cld(n, cols)
    xs = [((i - 1) % cols) + 1 for i in 1:n]
    ys = [((i - 1) ÷ cols) + 1 for i in 1:n]

    # --- Surface Z(row, col) = u ---
    Z = fill(NaN, rows, cols)
    for i in 1:n
        Z[ys[i], xs[i]] = zscale * u[i]
    end
    xgrid = collect(1:cols)
    ygrid = collect(1:rows)

    # 4 camera views: (azimuth, elevation)
    views = [
        (-pi/4,  pi/6),
        ( pi/4,  pi/6),
        (3pi/4,  pi/6),
        (-3pi/4, pi/6),
    ]

    fig = GLMakie.Figure(size=(1400, 1000))

    for k in 1:4
        r = (k <= 2) ? 1 : 2
        c = (k % 2 == 1) ? 1 : 2

        ax = GLMakie.Axis3(fig[r, c];
            title = "View $k",
            xlabel = "x", ylabel = "y", zlabel = "z = u",
            aspect = :data
        )

        Makie.surface!(ax, xgrid, ygrid, permutedims(Z))
        #Makie.scatter!(ax, xs, ys, zscale .* u; markersize=markersize)

        if length(path) > 1
            path_pts = GLMakie.Point3f.([xs[i] for i in path],
                                        [ys[i] for i in path],
                                        zscale .* [u[i] for i in path])
            Makie.lines!(ax, path_pts; linewidth=4, color=:red)
        end

        Makie.xlims!(ax, 0.5, cols + 0.5)
        Makie.ylims!(ax, 0.5, rows + 0.5)

        ax.azimuth[]   = views[k][1]
        ax.elevation[] = views[k][2]
    end

    filename = "$(filename)_3D_4views_cols_$(cols).png"
    if filename !== nothing
        Makie.save(filename, fig)
    end
    

    return fig, path
end