using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random

# Assuming square Ω!
# We also assume same numbering of nodes as in week2 
noelms = 3
function xy(noelms)
    Δ = 1 / noelms
    VX, VY = [], []

    for i ∈ 0:noelms, j ∈ 0:noelms
        push!(VX, i * Δ)
        push!(VY, 1 - j * Δ)
    end 

    return VX, VY
end 

VX, VY = xy(noelms)
println("VX = $VX")
println("VY = $VY")

# Defining the walls Γ
walls_predefined = [
    # outer wall 
    (1, 2), (2, 3), (3, 4), 
    (1, 5), (5, 9), (9, 13), 
    (4, 8), (8, 12), (12, 16), 
    (13, 14), (14, 15), (15, 16), 
    # maze part 
    (2, 6), (7, 11), 
    (10, 11), (10, 14)
]

function possible_walls(noelms, node_idx)
    A = zeros(Int64, noelms + 1, noelms + 1)
    for i in 1:(noelms+1)^2 
        A[i] = i 
    end 
    
    # padding of zeros around A 
    B = zeros(Int64, noelms + 1 + 2, noelms + 1 + 2) 
    B[2:end-1, 2:end-1] = A 

    r = findall(x -> x == node_idx, B)
    i, j = r[1][1], r[1][2]

    possible_node_idxs = [
        B[i - 1, j], 
        B[i, j + 1], 
        B[i + 1, j], 
        B[i, j - 1]
    ]
    legal_node_idxs = possible_node_idxs[findall(x -> x != 0, possible_node_idxs)]

    return legal_node_idxs
end

return possible_walls

function generate_walls(noelms, path_squares)
    walls = []
    # outer walls 

    # inner walls 

    return walls 
end

walls = walls_predefined
println("walls = $walls")

# Plots 
scatter(VX, VY, legend=nothing, color="blue")
for (i, j) ∈ walls
    plot!([VX[i], VX[j]], [VY[i], VY[j]], color="red")
end 

savefig("maze.png")