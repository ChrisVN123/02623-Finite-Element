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

println("EXERCISE 2.1")
println("a) CASE 1:")
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
println("VX = $VX")
println("VY = $VY")

println("b) CASE 1:")
display(conelmtab(4, 3))

println("\nEXERCISE 2.2")
println("CASE 2.2a:")
n = 4 # or 9? 
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
d, abc = basfun(n, VX, VY, EToV)

println("Δ = $d")
println("abc = ")
display(abc)

x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

println("CASE 2.2b:")
for k = 1:3 
    println("k = $k")
    n_vec = outernormal(n, k, VX, VY, EToV) # n = 4 
    display(n_vec)
end 

println("\nEXERCISE 2.3")
println("CASE 1:")
x0, y0 = 0, 0 
L1, L2 = 1, 1
noelms1, noelms2 = 4, 3 
qt = zeros(noelms1 * noelms2 * 2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)
display(A)

function spdiags(A)
    """
    A is square 
    Should work!
    """
    M = size(A)[1]
    B = zeros(M, M)
    d = []
    
    count = 1
    for k ∈ [-M:0; 1:M]
        kth_diag = diag(A, k)
        if !iszero(kth_diag)
            l = length(kth_diag)
            if k <= 0 
                B[1:l, count] = kth_diag 
            else 
                B[M-l+1:end, count] = kth_diag
            end 
            push!(d, k)
            count += 1
        end 
    end 

    return B[:, 1:count-1], d 
end

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)


println("CASE 2:")
q(x, y) = -6 * x + 2 * y - 2
x0, y0 = -2.5, -4.8
L1, L2 = 7.6, 5.9
noelms1, noelms2 = 4, 3

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
qt = q.(VX, VY)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)

println("\nEXERCISE 2.4")
println("CASE 1:")
x0, y0 = 0, 0 
L1, L2 = 1, 1
noelms1, noelms2 = 4, 3
qt = zeros(noelms1 * noelms2 * 2)
f = ones(noelms1 * noelms2 * 2)
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)
display(EToV)

A, b = assembly(VX, VY, EToV, lam1, lam2, qt)
bnodes = calculate_bnodes(noelms1, noelms2)
A, b = dirbc(bnodes, f, A, b)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)

println("Array(A)[1:13, 1:13] = ")
display(Array(A)[1:13, 1:13])

println("CASE 2:")
x0, y0 = -2.5, -4.8 
L1, L2 = 7.6, 5.9 
noelms1, noelms2 = 4, 3
q(x, y) = -6 * x + 2 * y - 2 
f24_2(x, y) = x^3 - x^2 * y + y^2 - 1
lam1, lam2 = 1, 1

VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
EToV = conelmtab(noelms1, noelms2)

A, b = assembly(VX, VY, EToV, lam1, lam2, q.(VX, VY))
bnodes = calculate_bnodes(noelms1, noelms2)
A, b = dirbc(bnodes, f24_2.(VX, VY), A, b)

B, d = spdiags(A)
println("B = ")
display(B)
println("d = ")
display(d)

println("b = ")
display(b)