using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random

######## 2.1 ########
function xy(x0, y0, L1, L2, noelms1, noelms2)
    Δx = L1 / noelms1 
    Δy = L2 / noelms2

    VX = Vector{Float64}(undef, (noelms1+1)*(noelms2+1))
    VY = Vector{Float64}(undef, (noelms1+1)*(noelms2+1))
    k=0
    for i ∈ 0:noelms1, j ∈ 0:noelms2
        k += 1
        VX[k] = x0 + i * Δx
        VY[k] = y0 + L2 - j * Δy
    end 

    return VX, VY
end

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

######## 2.2 ########
function basfun(n, VX, VY, EToV)
    abc = Matrix{Float64}(undef, 3, 3)
    
    i, j, k = EToV[n, :]
    x1, y1 = VX[i], VY[i]
    x2, y2 = VX[j], VY[j]
    x3, y3 = VX[k], VY[k]

    Δ = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

    # (i,j,k) = (1, 2, 3)
    abc[1, 1] = x2*y3 - x3*y2   # ã₁
    abc[1, 2] = y2 - y3         # b̃₁
    abc[1, 3] = x3 - x2         # c̃₁

    # (i,j,k) = (2, 3, 1)
    abc[2, 1] = x3*y1 - x1*y3   # ã₂
    abc[2, 2] = y3 - y1         # b̃₂
    abc[2, 3] = x1 - x3         # c̃₂

    # (i,j,k) = (3, 1, 2)
    abc[3, 1] = x1*y2 - x2*y1   # ã₃ 
    abc[3, 2] = y1 - y2         # b̃₃
    abc[3, 3] = x2 - x1         # c̃₃
    return Δ, abc
end

######## 2.3 ########
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

######## 2.4 ########
function calculate_bnodes(noelms1, noelms2)
    bnodes = []
    for i ∈ 0:noelms1, j ∈ 0:noelms2
        if (i == 0 || i == noelms1 || j == 0 || j == noelms2)
            # push!(bnodes, (i, j, j+1 + i * (noelms1 + 1)))
            push!(bnodes, j+1 + i * (noelms2 + 1))
        end 
    end 
    return bnodes
end

function dirbc(bnodes, f, A, b)
    """
    Algorithm 6 
    It seems that Γ₂ = Γ 

    bnodes = global node numbers of the nodes on Γ₂ 
    """

    M = size(A)[1]
    for i ∈ bnodes 
        A[i, i] = 1 
        b[i] = f[i] 
        for j ∈ 1:M
            if j != i 
                A[i, j] = 0 
                if !(j ∈ bnodes)
                    b[j] += -A[j, i] * f[i]
                    A[j, i] = 0
                end 
            end 
        end 
    end 

    return A, b 
end 

######## 2.6 ########
function local_node(col_idx)
    a = col_idx[1][1]
    b = col_idx[2][1]
    if sort([1,2]) == sort([a,b]) # true
        return 1
    elseif sort([2,3]) == sort([a,b])   # true
        return 2
    else return 3
    end
end

function ConstructBeds(noelms1, noelms2, EToV)
    bnodes = calculate_bnodes(noelms1, noelms2)
    E1 = Matrix{Int}(undef, 2*noelms1+2*noelms2, 2)
    k = 1
    noelms = noelms1 * noelms2 * 2
    noelms_to_exclusion = noelms2 * 2 - 1
    skip_vals = [noelms_to_exclusion, noelms - noelms_to_exclusion + 1] #check for other values
    for i in 1:size(EToV)[1]
        if i in skip_vals
            continue
        end
        intsct = intersect(Set(bnodes),Set(EToV[i,:]))
        l = length(intsct)
        global_nodes = []
        col_idx = []
        if l > 1
            for j in intsct
                push!(global_nodes, j)
                push!(col_idx, findall(x -> x == j, EToV[i,:]))
                
            end
            if l ==2
                idx = local_node(col_idx)
                E1[k,:] = [i, idx]
                k += 1
            else #l == 3
                E1[k,:] = [i, 3]
                k += 1
                E1[k,:] = [i, 1]
                k += 1
                
            end
        end
    end 
    return E1
end

######## 2.6 ########
function neubc(VX, VY, EToV, beds, q, b)
    """
    q is a function 
    """
    E1 = size(beds)[1]
    for p ∈ 1:E1 
        n, r = beds[p, 1], beds[p, 2]

        if r == 1 
            s = 2 
        elseif r == 2 
            s = 3 
        elseif r == 3 
            s = 1 
        end 
                
        i, j = EToV[n, r], EToV[n, s]
        xi, yi = VX[i], VY[i]
        xj, yj = VX[j], VY[j]

        xc, yc = (xi + xj) / 2, (yi + yj) / 2 
        q_midpoint = q(xc, yc)

        q_t = q_midpoint / 2 * ((xj - xi)^2 + (yj - yi)^2)^0.5 # since they are equal by (2.41)
        b[i] += -q_t
        b[j] += -q_t 
    end 

    return b 
end

######## 2.7 ########
function boundary_edges(noelms1, noelms2, beds)
    Γ1 = Matrix{Int}(undef, noelms1+noelms2, 2)
    Γ2 = Matrix{Int}(undef, noelms1+noelms2, 2)

    k1 = 1
    k2 = 1
    for (i, v) in enumerate(beds[:, 1])
        if v % 2 == 0 
            Γ1[k1, 1] = v
            Γ1[k1, 2] = beds[i, 2]
            k1 += 1
        elseif v % 2 == 1 
            Γ2[k2, 1] = v
            Γ2[k2, 2] = beds[i, 2]
            k2 += 1
        end 
    end 
    
    return Γ1, Γ2
end

function dirichlet_bound(noelms1, noelms2)   
    bnodes_dirbc = []
    for i in 0:(noelms1)
        push!(bnodes_dirbc, 1 + i * (noelms2+1))
    end 
    bnodes_dirbc = [bnodes_dirbc; collect(last(bnodes_dirbc)+1:last(bnodes_dirbc)+noelms2)]
    
    return bnodes_dirbc
end

function compute_error(VX, VY, û, u)
    """ 
    u is a function! 
    Estimating the error between the solution and the approximation at the nodes
    """
    return maximum(
        abs.(û - u.(VX, VY))
    )
end

######## 2.8 ########
function Driver28b(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, q̃)
    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV = conelmtab(noelms1, noelms2)

    u_x(x, y) = 0
    u_y(x, y) = 0
    q_l(x, y) = u_x(x, y)
    q_b(x, y) = u_y(x, y)

    # GLOBAL ASSEMBLY 
    A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
    # bnodes = calculate_bnodes(noelms1, noelms2)
    beds = ConstructBeds(noelms1, noelms2, EToV)
    Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

    # IMPOSING NEUMANN 
    b_l = neubc(VX, VY, EToV, Γ1[1:noelms2, :], q_l, b)
    b_b = neubc(VX, VY, EToV, Γ1[noelms1:end, :], q_b, b_l) 

    # IMPOSING DIRICHLET 
    bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
    A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b_b)

    # SOLUTION 
    û = A \ b

    return VX, VY, EToV, û
end

function Driver28c(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, q̃)
    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV = conelmtab(noelms1, noelms2)

    # GLOBAL ASSEMBLY 
    A, b = assembly(VX, VY, EToV, lam1, lam2, q̃.(VX, VY))
    # bnodes = calculate_bnodes(noelms1, noelms2)
    beds = ConstructBeds(noelms1, noelms2, EToV)
    Γ1, Γ2 = boundary_edges(noelms1, noelms2, beds)

    # IMPOSING DIRICHLET 
    bnodes_dirbc = dirichlet_bound(noelms1, noelms2)
    A, b = dirbc(bnodes_dirbc, f.(VX, VY), A, b)

    # SOLUTION
    û = A \ b

    return VX, VY, EToV, û
end