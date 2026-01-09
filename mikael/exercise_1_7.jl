using Plots 
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random
using Polynomials
# functions

#from 1.2
function K(h_list, r, s, i)
    if (r == 1 && s == 1) || (r == 2 && s == 2)
        return 1 / h_list[i] + h_list[i] / 3
    else 
        return -1 / h_list[i] + h_list[i] / 6
    end 
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

function compute_error_decrease(fun, VX, EToV)
    M = length(VX)
    err_list = []
    num = 1_000 # must be even! 
    for (x_i, x_ip1) in zip(VX[1:end-1], VX[2:end])
        x_i_new = x_i + (x_ip1 - x_i) / 2 
        x_arr = LinRange(x_i, x_ip1, num)
        h = (x_ip1 - x_i) / num 

        u_I_h = u_I.(x_arr, [[x_i, x_ip1]])

        half = div(num, 2) 

        u_I_h_refined_1 = u_I.(x_arr[1:half-1], [[x_i, x_i_new]])
        u_I_h_refined_2 = u_I.(x_arr[half:end], [[x_i_new, x_ip1]])
        u_I_h_refined = [u_I_h_refined_1; u_I_h_refined_2]

        # approximate the ||.||_2 by the sum 
        err = sqrt(sum(
            (u_I_h - u_I_h_refined).^2 .* h 
        )) 
        # use analytical solution! 

        err_list = [err_list; err]
    end 

    return err_list
end

# b) 
print("b)\n")
# Still NOT using EToVcoarse 
function refine_marked(EToVcoarse, xcoarse, idxMarked)
    EToVfine = EToVcoarse
    xfine = xcoarse
    
    for (i, idx) in enumerate(idxMarked)
        if idx == 1 
            x_i, x_ip1 = xcoarse[i], xcoarse[i+1]
            dist = (x_ip1 - x_i) / 2
            x_new = x_i + dist 
            
            insert!(xfine, i+1, x_new)
            insert!(idxMarked, i+1, 0) # this is to take into account that xfine grows which makes xcoarse grow as well.
        end
    end 

    return EToVfine, xfine
end 

function calc()
    idxMarked = [1, 1]
    xfine = copy(xcoarse_initial)
    while sum(idxMarked) != 0
        idxMarked = Int.(
            compute_error_decrease("", xfine, "") .> 10^-4
        )
        _, xfine = refine_marked("", xfine, idxMarked)
    end 

    return xfine
end 

function u_hat(x, u_coeffs_hat, VX)
    res = 0 
    for (i, u_coeff_hat) in enumerate(u_coeffs_hat)
        res += u_coeff_hat * N(x, i, VX)
    end 
    return res 
end 

"""
We consider the BVP 
    u''(x) - u(x) = f(x), 0 <= x <= 1 
    u(0) = c, u(1) = d 
"""

function u(x)
    return exp(-800*(x - 0.4)^2) + 0.25*exp(-40*(x - 0.8)^2)
end 

# a) 
function f(x)
    return exp(-32*(5*x - 2)^2)*(2560000*x^2 - 2048000*x + 407999) + exp(-(8*(5*x - 4)^2)/5)*(-81/4 + 64*(5*x - 4)^2)
end 

xc = [0, 0.5, 1]
xf = [0, 0.25, 0.5, 1]
uhc = [10,20,10]
uhf = [9,14,18,10]
EToVc = [[1,2],[2,3]]               # not used! 
EToVf = [[1,2],[2,3],[3,4],[4,5]]   # not used!

# b) 
print("b)\n")
function create_mapping_coarse_fine(xc, xf)
    res = []
    for x in xc, (j, y) in enumerate(xf)
        if x == y 
            res = [res; j]
        end 
    end 
    return res 
end 

function error_estimate_working(xc, xf, uhc, uhf)
    Old2New = create_mapping_coarse_fine(xc, xf)

    err_arr = []
    for (i, (x_i, x_ip1)) in enumerate(zip(xc[1:end-1], xc[2:end]))
        x_arr_arr = []
        uhf_arr_arr = []
        
        idx_lower, idx_upper = Old2New[i], Old2New[i+1]
        x_lower = x_i 
        for j in idx_lower:(idx_upper-1) 
            x_upper = xf[j+1]
            x_arr = LinRange(x_lower, x_upper, 500)
            uhf_arr = u_hat.(x_arr, [uhf], [xf])

            x_arr_arr = [x_arr_arr; x_arr]
            uhf_arr_arr = [uhf_arr_arr; uhf_arr]

            x_lower = x_upper
        end 

        uhc_arr_arr = u_hat.(x_arr_arr, [uhc], [xc])
        h = (x_ip1 - x_i) / length(x_arr_arr) # approximation! 
        err = sqrt(sum(
            (uhc_arr_arr .- uhf_arr_arr).^2 * h 
        ))
        err_arr = [err_arr; err]
    end

    return err_arr
end

# xc = [0, 0.5, 1.0]
# xf = [0, 0.25, 0.5, 0.75, 1.0]
# uhc = [0, 10, 5]
# uhf = [0, 3, 10, 11, 5]
# error_estimate_working(xc, xf, uhc, uhf)

#c)
print("c)\n")
function refine_marked(EToVcoarse, xcoarse, idxMarked)
    EToVfine = EToVcoarse
    xfine = copy(xcoarse)
    
    j = 0
    for (i, idx) in enumerate(idxMarked)
        if idx == 1 
            x_i, x_ip1 = xcoarse[i], xcoarse[i+1]
            dist = (x_ip1 - x_i) / 2
            x_new = x_i + dist 
            
            insert!(xfine, i+1+j, x_new)
            j += 1
        end
    end 

    return EToVfine, xfine
end 

# d) 
# function create_b(f_list, x)
#     M = length(f_list)
#     h_list = x[2:end] - x[1:end-1]
#     b = zeros(M)
    
#     for i = 2:(M-1)
#         b[i] = f_list[i-1] * h_list[i-1] / 6 + f_list[i] * (h_list[i-1] + h_list[i]) / 3 + f_list[i+1] * h_list[i] / 6
#     end     
    
#     return b 
# end 

function create_b(M, x)
    h_list = x[2:end] - x[1:end-1]

    d = 150

    b = zeros(M)
    for i = 2:(M-1)
        x_im1, x_ip1 = x[i-1], x[i+1]
        x_arr = LinRange(x_im1, x_ip1, d)
        h = (x_ip1 - x_im1) / d
        
        y_arr = f.(x_arr) .* N.(x_arr, [i], [x])
        b[i] = sum(y_arr) * h
    end 

    return b 
end 

function BVP1Drhs(L, c, d, x, func)
    M = length(x)
    h_list = x[2:end] - x[1:end-1]

    # Algorithm 1
    rows = Int.(zeros(4 * M))
    cols = Int.(zeros(4 * M))
    vals = zeros(4 * M)
    # f_list = f.(x)
    # b = -create_b(f_list, x)
    b = -create_b(M, x)
    
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
end 

print("d)\n")
print("e) and f)\n")
function calc2(L, c, d, xc)
    M = length(xc)

    Δerr_i = 10^(-4)
    idxMarked = zeros(M-1)
    idxMarked[M-1] = 1
    idxMarked = Int.(idxMarked)
    max_N = 1_000
    k = 0
    xf = 0 

    x_arr = LinRange(0, 1, 2_000)
    h_arr = []
    err_arr = []

    while sum(idxMarked) != 0 && k <= max_N
        print("$k, $idxMarked\n")
        _, xf = refine_marked("", xc, idxMarked)

        _, _, uhc = BVP1Drhs(L, c, d, xc, f)
        _, _, uhf = BVP1Drhs(L, c, d, xf, f)
        error_est = error_estimate_working(xc, xf, uhc, uhf) 
        print("$(length(xc)), $(length(xf))\n")

        ##############################################################
        ################# ONLY USED FOR CONVERGENCE ##################
        # error bounds a la 1.7 
        err = maximum(abs.(
            u.(x_arr) - u_hat.(x_arr, [uhf], [xf])
        ))
        h = maximum(xf[2:end] - xf[1:end-1])
        h_arr = [h_arr; h]
        err_arr = [err_arr; err]
        ##############################################################

        idxMarked = Int.(error_est .> Δerr_i)
        xc = copy(xf)
        k += 1
    end 
    
    return k, xf, h_arr, err_arr
end 

xc = [0, 0.5, 1.0]
# xc = collect(LinRange(0, 1, 50))
L = 1
c = exp(-128) + 1 / 4 * exp(-128 / 5)
d = exp(-288) + exp(-8 / 5)/4
print("c, d = $c, $d\n")

k, xf, h_arr, err_arr = calc2(L, c, d, xc)
_, _, u_fine = BVP1Drhs(L, c, d, xf, f)

x_arr = LinRange(0, 1, 1_000)
# x_arr = []
# for (x_i, x_ip1) in zip(xf[1:end-1], xf[2:end])
#     # print("$x_i, $x_ip1\n")
#     # x_arr = [x_arr; LinRange(x_i, x_ip1, 1_000)]
#     append!(x_arr, LinRange(x_i, x_ip1, 3))
# end 
u_arr = u.(x_arr)
u_fine_arr = u_hat.(x_arr, [u_fine], [xf])

plot(x_arr, u_fine_arr, label="u_hat", title="u vs u_hat")
plot!(x_arr, u_arr, label="u")
xlims!(0, 1)
# ylims!(-12, 1)

savefig("exercise_1_7_e_u_hat_xf.png")

print("length(xf) = $(length(xf))\n")
print("k iterations: $k\n")

# print("\n$h_arr\n$err_arr\n\n")
print("fit: $(fit(log.(h_arr), log.(err_arr), 1))\n")

histogram(xf)
savefig("exercise_1_7_histogram.png")