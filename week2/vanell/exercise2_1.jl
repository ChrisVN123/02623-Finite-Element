using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

function xy(x0, y0, L1, L2, noelms1, noelms2)
    Δx = L1 / noelms1 
    Δy = L2 / noelms2

    VX = Vector{Float64}(undef, (noelms1+1)*(noelms2+1))
    VY = Vector{Float64}(undef, (noelms1+1)*(noelms2+1))
    k=0
    for i in 0:noelms1, j in 0:noelms2
        k += 1
        VX[k] = x0 + i * Δx
        VY[k] = y0 + L2 - j * Δy
    end 

    return VX, VY
end

x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
println("Exercise 2.1 test case 1")
println("VX = $VX")
println("VY = $VY")
function conelmtab(nx::Int, ny::Int)
    @assert nx ≥ 2 && ny ≥ 2 
    nx += 1
    ny += 1
    num_tris = 2 * (nx - 1) * (ny - 1)
    etoV = Matrix{Int}(undef, num_tris, 3)

    e = 1
    tl = 1 

    for col in 1:(nx - 1)            
        tl = 1 + (col - 1) * ny   
        for row in 1:(ny - 1)          
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

println("Exercise 2.1 b:")
display(conelmtab(4, 3))
