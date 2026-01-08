using Plots 
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials

function K(h_list, r, s, i)
    if (r == 1 && s == 1) || (r == 2 && s == 2)
        return 1 / h_list[i] + h_list[i] / 3
    else 
        return -1 / h_list[i] + h_list[i] / 6
    end 
end 

# a) 
print("a)\n")
function BVP1D_a(L, c, d, x)
    M = length(x)
    h_list = x[2:end] - x[1:end-1]

    # Algorithm 1
    rows = Int.(zeros(4 * M))
    cols = Int.(zeros(4 * M))
    vals = zeros(4 * M)
    b = zeros(M)
    
    count = 0
    for i in 1:(M-1)
        rows[count .+ (1:4)] = [i i   i+1 i+1]
        cols[count .+ (1:4)] = [i i+1 i+1 i  ]
        vals[count .+ (1:4)] = [
            K(h_list, 1, 1, i) 
            K(h_list, 1, 2, i)
            K(h_list, 2, 2, i)
            K(h_list, 2, 1, i)
        ]

        count += 4
    end 
    A = sparse(rows[1:count], cols[1:count], vals[1:count])

    # Algorithm 2
    ## NOTE: check if A[1, 2] is efficient 
    b[1] = c
    b[2] = b[2] - A[1, 2] * c 
    A[1, 1] = 1
    A[1, 2] = 0 
    A[2, 1] = 0 # modified 
    b[M] = d 
    b[M - 1] = b[M - 1] - A[M - 1, M] * d
    A[M, M] = 1 
    A[M - 1, M] = 0
    A[M, M - 1] = 0 # modified

    # Version 1
    return A, b, A \ b

    # Version 2 
    # return A, b, cholesky(A) \ b
end 

L = 2 
c = 1 
d = exp(2)
x = [0.0, 0.2, 0.4, 0.6, 0.7, 0.9, 1.4, 1.5, 1.8, 1.9, 2.0]

BVP1D_a(L, c, d, x)
print("BVP1D_a timing:\n")
@time BVP1D_a(L, c, d, x)
A, b, u_coeffs_hat = BVP1D_a(L, c, d, x)

x_arr = LinRange(0, L, 1_000)
u_correct = exp.(x_arr)
plot(x_arr, u_correct, label = "u")
plot!(x, u_coeffs_hat, label = "hat(u)")
title!("Exercise 1.2 a) jl")
xlims!(0, 2)
ylims!(1, 7)
savefig("exercise_1_2_a_jl.png")

# b) 
print("b)\n")
function BVP1D_b(L, c, d, M)
    h = L / M 
    k11 = 1 / h + h / 3 
    k12 = -1 / h + h / 6
    
    # Algorithm 1
    rows = Int.(zeros(4 * M))
    cols = Int.(zeros(4 * M))
    vals = zeros(4 * M)
    b = zeros(M)
    
    count = 0
    for i in 1:(M-1)
        rows[count .+ (1:4)] = [i i   i+1 i+1]
        cols[count .+ (1:4)] = [i i+1 i+1 i  ]
        vals[count .+ (1:4)] = [k11 k12 k11 k12]

        count += 4
    end 
    A = sparse(rows[1:count], cols[1:count], vals[1:count])

    # Algorithm 2
    ## NOTE: check if A[1, 2] is efficient 
    b[1] = c
    b[2] = b[2] - A[1, 2] * c 
    A[1, 1] = 1
    A[1, 2] = 0 
    A[2, 1] = 0 # modified 
    b[M] = d 
    b[M - 1] = b[M - 1] - A[M - 1, M] * d
    A[M, M] = 1 
    A[M - 1, M] = 0
    A[M, M - 1] = 0 # modified

    return A, b, A \ b
end 

L = 2 
c = 1 
d = exp(2)
M = 11 

BVP1D_b(L, c, d, M)
print("BVP1D_b timing:\n")
@time BVP1D_b(L, c, d, M)
A, b, u_coeffs_hat = BVP1D_b(L, c, d, M)
print(u_coeffs_hat)

VX = LinRange(0, L, M)

plot(x_arr, u_correct, label = "u")
plot!(VX, u_coeffs_hat, label = "hat(u)")
title!("Exercise 1.2 b) jl")
xlims!(0, 2)
ylims!(1, 7)
savefig("exercise_1_2_b_jl.png")

# c) 
print("c)\n")
## Looks good right? 

# d) 
print("d)\n")
function u(x)
    return exp(x)
end 

function N(x, i, VX)
    M = length(VX)
    if i == 1 
        x1, x2 = VX[1], VX[2]
        if x1 <= x && x <= x2 
            1 - (x - x1) / (x2 - x1)
        else 0 
        end 
    elseif i == M 
        x_Mm1, x_M = VX[M-1], VX[M]
        if x_Mm1 <= x && x <= x_M
            (x - x_Mm1) / (x_M - x_Mm1)
        else 0
        end 
    elseif 2 <= i && i <= (M - 1)
        x_im1, x_i, x_ip1 = VX[i-1], VX[i], VX[i+1]
        if x_im1 <= x && x <= x_i 
            (x - x_im1) / (x_i - x_im1)
        elseif x_i <= x && x <= x_ip1
            1 - (x - x_i) / (x_ip1 - x_i) 
        else 0 
        end 
    end 
end

function u_hat(x, u_coeffs_hat, VX)
    res = 0 
    for (i, u_coeff_hat) in enumerate(u_coeffs_hat)
        res += u_coeff_hat * N(x, i, VX)
    end 
    return res 
end 

function error_calculation(M)
    VX = LinRange(0, L, M)
    A, b, u_coeffs_hat = BVP1D_b(L, c, d, M)
    x_arr = LinRange(0, L, 2_000) # can be changed 

    error = maximum(abs.(
        u_hat.(x_arr, [u_coeffs_hat], [VX]) .- u.(x_arr)
    ))

    return error 
end 

M_arr = 2_000:100:5_000
h_arr = L ./ M_arr 
error_arr = error_calculation.(M_arr)

plot(h_arr, error_arr)
savefig("exercise_1_2_d_jl.png")

# since error_arr \approx C h_arr^2 then 
# log(error_arr) \approx 2 log(h_arr) + log(C)
# so (polynomial)fit!

print("polynomial fit: \n$(fit(log.(h_arr), log.(error_arr), 1))\n")