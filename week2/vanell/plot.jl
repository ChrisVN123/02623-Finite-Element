using Pkg

# Pkg.add("GLMakie")
# Pkg.add("GeometryBasics")

using GLMakie
using GeometryBasics

include("exercise2_1.jl")
include("exercise2_2.jl")
include("exercise2_3.jl")
include("exercise2_4.jl")


function plot_fem_solution3d(VX, VY, EToV, u;
    colormap=:viridis,
    show_wireframe::Bool=true,
    wire_color=(:black, 0.25),
    figure_size=(900, 700),
    title::AbstractString="FEM solution",
    colorrange=nothing,
    filename::AbstractString="plottingFEM.png"
)
    @assert length(VX) == length(VY)
    @assert length(u)  == length(VX)[1]
    @assert ndims(EToV) == 2

    x = Float32.(VX)
    y = Float32.(VY)
    z = Float32.(u)

    nen = size(EToV, 2)
    @assert nen == 3 || nen == 4

    positions = GeometryBasics.Point3f.(x, y, z)
    

    faces = GeometryBasics.TriangleFace{Int}[]
    if nen == 3
        for e in eachrow(EToV)
            a, b, c = e
            push!(faces, GeometryBasics.TriangleFace(a, b, c))
        end
    else
        for e in eachrow(EToV)
            a, b, c, d = e
            push!(faces, GeometryBasics.TriangleFace(a, b, c))
            push!(faces, GeometryBasics.TriangleFace(a, c, d))
        end
    end

    msh = GeometryBasics.Mesh(positions, faces)

    fig = Figure(size=figure_size)
    ax = Axis3(fig[1, 1], title=title, xlabel="x", ylabel="y", zlabel="u")

    cr = colorrange === nothing ? (minimum(z), maximum(z)) : colorrange
    plt = mesh!(ax, msh; color=z, colormap=colormap, colorrange=cr)

    if show_wireframe
        GLMakie.wireframe!(ax, msh; color=wire_color)
    end

    Colorbar(fig[1, 2], plt, label="u")

    # Saves in the current working directory:
    save(filename, fig)

    return fig
end


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


# x0, y0 = -2.5, -4.8 
# L1, L2 = 7.6, 5.9 
# noelms1, noelms2 = 4, 3
# q(x, y) = -6 * x + 2 * y - 2 
# f24_2(x, y) = x^3 - x^2 * y + y^2 - 1
# lam1, lam2 = 1, 1

# VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
# print(length(VX))

# EToV = conelmtab(noelms1, noelms2)
# print(size(EToV, 2))
# A, b = assembly(VX, VY, EToV, lam1, lam2, q.(VX, VY))
# bnodes = calculate_bnodes(noelms1, noelms2)
# A, b = dirbc(bnodes, f24_2.(VX, VY), A, b)

# u = A \ b
# fig = plot_fem_solution3d(VX, VY, EToV, u)  # writes ./plottingFEM.png
# save(joinpath(pwd(), "plottingFEM.png"), fig)