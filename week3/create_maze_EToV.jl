using Random

function maze_EToV(rows::Int, cols::Int; seed=nothing, braid::Float64=0.0, startcell::Tuple{Int,Int}=(1,1))


    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

    idx(r, c) = (r - 1) * cols + c
    inbounds(r, c) = (1 ≤ r ≤ rows) && (1 ≤ c ≤ cols)

    # Carved connectivity (adjacency sets) for degree checks / braiding
    adj = [Int[] for _ in 1:(rows * cols)]

    visited = falses(rows, cols)
    stack = Tuple{Int,Int}[]

    sr, sc = startcell

    push!(stack, (sr, sc))
    visited[sr, sc] = true

    # Helper: list unvisited neighbors (4-neighborhood)
    function unvisited_neighbors(r, c)
        nbrs = Tuple{Int,Int}[]
        for (dr, dc) in ((-1,0), (1,0), (0,-1), (0,1))
            rr, cc = r + dr, c + dc
            if inbounds(rr, cc) && !visited[rr, cc]
                push!(nbrs, (rr, cc))
            end
        end
        return nbrs
    end

    edges = Tuple{Int,Int}[]
    function carve!(r1, c1, r2, c2)
        a = idx(r1, c1)
        b = idx(r2, c2)
        if a > b
            a, b = b, a
        end
        push!(edges, (a, b))
        push!(adj[idx(r1,c1)], idx(r2,c2))
        push!(adj[idx(r2,c2)], idx(r1,c1))
        return nothing
    end

    #Randomized DFS to backtrack ap path
    while !isempty(stack)
        r, c = stack[end]
        nbrs = unvisited_neighbors(r, c)
        if isempty(nbrs)
            pop!(stack)
        else
            (rr, cc) = rand(rng, nbrs)
            carve!(r, c, rr, cc)
            visited[rr, cc] = true
            push!(stack, (rr, cc))
        end
    end

    # Optional braiding for add extra edges at dead ends ---
    if braid > 0.0
        for r in 1:rows, c in 1:cols
            v = idx(r, c)
            if length(adj[v]) == 1 && rand(rng) < braid
                # candidate neighbors not already connected
                candidates = Int[]
                for (dr, dc) in ((-1,0), (1,0), (0,-1), (0,1))
                    rr, cc = r + dr, c + dc
                    if inbounds(rr, cc)
                        u = idx(rr, cc)
                        if !(u in adj[v])
                            push!(candidates, u)
                        end
                    end
                end
                if !isempty(candidates)
                    u = rand(rng, candidates)
                    a, b = min(v, u), max(v, u)
                    push!(edges, (a, b))
                    push!(adj[v], u)
                    push!(adj[u], v)
                end
            end
        end
    end

    #remove any accidental duplicates (can occur with braiding choices)
    sort!(edges)
    unique!(edges)

    return edges
end
