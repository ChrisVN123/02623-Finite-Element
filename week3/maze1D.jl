using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random

include("plotpath_2d.jl")
include("plotpath_3d.jl")
include("EToV_50x50_from_image.jl")
include("plotpath.jl")
include("create_maze_EToV.jl")

###inputs
start, stop, elements = 25, 2476, 2500
cols = 50
c = 1
d = 0

function assemble_K(EToV, c, d, start, stop, n_elements)
    E = length(EToV)
    b = zeros(n_elements)

    I = Int[]
    J = Int[]
    V = Float64[]

    for e in 1:E
        i, j = EToV[e]

        push!(I, i); push!(J, i); push!(V,  1)
        push!(I, j); push!(J, j); push!(V,  1)
        push!(I, i); push!(J, j); push!(V, -1)
        push!(I, j); push!(J, i); push!(V, -1)
    end
    A =  sparse(I, J, V, n_elements, n_elements)


    # Algorithm 2
    b[start] = c
    b[stop] = d
    
    # Apply boundary conditions
    A[start, :] .= 0
    A[stop, :] .= 0
    A[start, start] = 1
    A[stop, stop] = 1
    
    #don't explicity write zeros
    A = dropzeros(A)    
    return A, b
end

# A,b = (assemble_K(EToV, c, d,  start, stop, elements))
# @time u = A \ b
# @time A \ b

# @time Array(A) \ b
# plotting_path2d(EToV, u, cols, start)
# plotting_path_3d(EToV, u, cols, start;  zscale = Float64(10))
# plotting_path_3d_4views(EToV, u, cols, start;  zscale = Float64(10))
# plotting_path(EToV, u, cols, start)

rows, cols = 5, 4
EToV = maze_EToV(rows, cols; seed=42, braid=0.05)   # braid optional

elements = rows * cols
start = 1
stop  = 20
c, d = 1, 0

A, b = assemble_K(EToV, c, d, start, stop, elements)

u = A \ b

plotting_path2d(EToV, u, cols, start)
plotting_path_3d(EToV, u, cols, start; zscale=10.0)
plotting_path_3d_4views(EToV, u, cols, start; zscale=10.0)
plotting_path(EToV, u, cols, start)
