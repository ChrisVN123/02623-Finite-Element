function conelmtab(nx::Int, ny::Int)
    @assert nx ≥ 2 && ny ≥ 2 
    nx += 1
    ny += 1
    num_tris = 2 * (nx - 1) * (ny - 1)
    etoV = Matrix{Int}(undef, num_tris, 3)

    e = 1
    tl = 1 

    for col ∈ 1:(nx - 1)            
        tl = 1 + (col - 1) * ny   
        for row ∈ 1:(ny - 1)          
            tr = tl + ny
            bl = tl + 1
            br = tl + ny + 1
            
            #tri1
            etoV[e, 1] = tr
            etoV[e, 2] = tl
            etoV[e, 3] = br
            e += 1

            #tri2
            etoV[e, 1] = bl
            etoV[e, 2] = br
            etoV[e, 3] = tl
            e += 1

            tl += 1  #move one node down within the column
        end
    end

    return etoV
end

function assembly(VX, VY, EToV, lam1, lam2, qt)
    """
    Algorithm 4 
    """
    
    N = size(EToV)[1]
    M = length(VX)
    
    A = zeros(Float64, M, M)
    b = zeros(Float64, M)  # vector is simpler/safer

    for n ∈ 1:N
        i, j, k = EToV[n, :]
        Δ, abc = basfun(n, VX, VY, EToV) 
        
        for r ∈ 1:3
            q̃ = (qt[i] + qt[j] + qt[k]) / 3 # Can be moved out
            q_r_n = abs(Δ)/3 * q̃
            ĩ = EToV[n, r]
            
            b[ĩ] += q_r_n
            for s ∈ 1:3
                j̃ = EToV[n, s]

                b_r, b_s = abc[r, 2], abc[s, 2]
                c_r, c_s = abc[r, 3], abc[s, 3]
                
                k_rs = 1 / (4 * abs(Δ)) * (lam1 * b_r * b_s + lam2 * c_r * c_s)
                A[ĩ, j̃] += k_rs
            end 
        end 
    end

    A = sparse(A)
    return A, b 
end 