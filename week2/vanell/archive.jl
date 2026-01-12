values = Matrix{NTuple{2,Float64}}(undef, n, 3)

for e in 1:n, j in 2:4
    vid = EToV[e, j]  # vertex id
    idx = findfirst(t -> t[end] == vid, coordinates)
    idx === nothing && error("Vertex id $vid not found in coordinates.") #idx of vertex
    row, col = Tuple(idx)
    x, y = coordinates[row, col][1:2]
    values[e, j-1] = (x, y)
end