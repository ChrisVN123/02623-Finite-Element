
function bfs(adj)
    V = length(adj)
    visited = falses(V)
    res = Int[]

    src = 0
    q = Vector{Int}()
    visited[src + 1] = true          # Julia is 1-indexed
    push!(q, src)

    while !isempty(q)
        curr = popfirst!(q)          # dequeue from front
        push!(res, curr)

        # visit all unvisited neighbors
        for x in adj[curr + 1]       # adjust index for adjacency list
            if !visited[x + 1]
                visited[x + 1] = true
                push!(q, x)
            end
        end
    end

    return res
end
