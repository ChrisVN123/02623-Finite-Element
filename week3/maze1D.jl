using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random



function EToV_radius(nodes::AbstractVector, r::Real)
    n = length(nodes)
    r2 = r^2
    edges = Tuple{Int,Int}[]

    for i in 1:n-1
        xi, yi = nodes[i]
        for j in i+1:n
            xj, yj = nodes[j]
            if (xi - xj)^2 + (yi - yj)^2 ≤ r2
                push!(edges, (i, j))
            end
        end
    end
    return edges
end

Nodes = [[0,0], [1,0], [2,0]]#,[3,0], [1,-1]]
EToV = EToV_radius(Nodes, 1.0)
display(EToV)  # [(1, 2), (2, 3), (2, 4)]


function assemble_K(nodes, EToV, k)
    N = length(nodes)
    E = length(EToV)

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

    return sparse(I, J, V, N, N)
end


k = [1, 1, 1]#,1]  # define these scalars
A = collect(assemble_K(Nodes, EToV, k))
println("See A matrix")
display(collect(A))

# A_solve = collect(A)[2:4,2:4]
# b_solve = [0 0 0]
# display(b_solve)
# display(A\b_solve)

print(det(A))