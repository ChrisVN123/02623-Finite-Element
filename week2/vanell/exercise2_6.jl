using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

include("exercise2_1.jl")
include("exercise2_2.jl")
include("exercise2_3.jl")
include("exercise2_4.jl")

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

# x0, y0 = -2.5, -4.8 
# L1, L2 = 7.6, 5.9 
# noelms1, noelms2 = 4, 3
# q(x, y) = -6 * x + 2 * y - 2 
# f24_2(x, y) = x^3 - x^2 * y + y^2 - 1
# lam1, lam2 = 1, 1

# VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
# EToV = conelmtab(noelms1, noelms2)
# A, b = assembly(VX, VY, EToV, lam1, lam2, q.(VX, VY))

# beds = ConstructBeds(VX, VY, EToV)
# println("beds = ")
# display(beds)

function neubc(VX, VY, EToV, beds, q, b)
    """
    q is a function 
    """
    # println("b_input $b")
    E1 = size(beds)[1]
    for p ∈ 1:E1 
        # println("p, $p")
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

        # q_midpoint = q[p] # discuss! 
        xc, yc = (xi + xj) / 2, (yi + yj) / 2 
        q_midpoint = q(xc, yc)

        q_t = q_midpoint / 2 * ((xj - xi)^2 + (yj - yi)^2)^0.5 # since they are equal by (2.41)
        b[i] += -q_t
        b[j] += -q_t 
    end 

    return b 
end

# println(neubc(VX,VY,EToV, cb, q, b))