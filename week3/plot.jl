

using Graphs
using CairoMakie




EToV = [[1,2], [2,3], [2,5], [5,4], [4,6], [6,7], [7, 8],[8, 9]]

# ---- YOU SET THIS ----
cols = 3   # number of nodes in each row (set this to what you have)
# ----------------------

# Build graph
n = maximum(vcat(EToV...))
g = SimpleGraph(n)
for e in EToV
    add_edge!(g, e[1], e[2])
end

rows = cld(n, cols)

# Grid coordinates for each node (1-indexed)
xs = [((i - 1) % cols) + 1 for i in 1:n]
ys = [rows - ((i - 1) ÷ cols) for i in 1:n]  # flip y so row 1 is at top

fig = Figure(size=(900, 700))
ax = Axis(fig[1, 1];
    aspect = DataAspect(),
    xticks = 1:cols,
    yticks = 1:rows,
    title = "Maze graph on grid"
)

# draw grid
for x in 1:cols
    vlines!(ax, x, color=:gray90)
end
for y in 1:rows
    hlines!(ax, y, color=:gray90)
end

# draw edges
for e in EToV
    u, v = e
    lines!(ax, [xs[u], xs[v]], [ys[u], ys[v]])
end

# draw nodes
scatter!(ax, xs, ys, markersize=12)

# labels (optional)
text!(ax, string.(1:n), position=Point2f.(xs .+ 0.1, ys .+ 0.1), fontsize=12)

xlims!(ax, 0.5, cols + 0.5)
ylims!(ax, 0.5, rows + 0.5)

save("maze_grid.png", fig)

