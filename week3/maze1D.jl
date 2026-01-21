using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random

include("plotpath_2d.jl")
include("plotpath_3d.jl")
include("EToV_50x50_from_image.jl")
include("plotpath.jl")

###inputs
start, stop, elements = 2, 8, 9
cols = 3
c = 1
d = 0

function assemble_K(EToV, c, d, start, stop, n_elements)
    E = length(EToV)
    b = zeros(n_elements)
    k = ones(E)

    I = Int[]
    J = Int[]
    V = Float64[]

    for e in 1:E
        i, j = EToV[e]
        ke = float(k[e])

        push!(I, i); push!(J, i); push!(V,  ke)
        push!(I, j); push!(J, j); push!(V,  ke)
        push!(I, i); push!(J, j); push!(V, -ke)
        push!(I, j); push!(J, i); push!(V, -ke)
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

#bfs_adj =  bfs(EToV)
#print(bfs_adj)

A,b = (assemble_K(EToV, c, d,  start, stop, elements))
@time u = A \ b

@time A \ b

@time Array(A) \ b
plotting_path2d(EToV, u, cols, start)
plotting_path_3d(EToV, u, cols, start;  zscale = Float64(10))
plotting_path_3d_4views(EToV, u, cols, start;  zscale = Float64(10))
plotting_path(EToV, u, cols, start)