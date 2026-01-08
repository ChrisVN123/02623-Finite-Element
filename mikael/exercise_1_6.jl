using Plots
using BenchmarkTools

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

function u(x)
    exp(-800 * (x - 0.4)^2) + 0.25 * exp(-40 * (x - 0.8)^2)
end

function u_I(x, VX)
    M = length(VX)
    u_i_list = u.(VX)
    res = 0
    for i in 1:M
        res += u_i_list[i] * N(x, i, VX)
    end 

    return res
end

# a) 
print("a)\n")
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

# c) 
print("c)\n")
xcoarse_initial = [0, 0.5, 1.0]
delta_err_i = 10^-4

function calc()
    idxMarked = [1, 1]
    xfine = copy(xcoarse_initial)
    count = 0 
    while sum(idxMarked) != 0 && count < 500
        idxMarked = Int.(
            compute_error_decrease("", xfine, "") .> delta_err_i
        )
        _, xfine = refine_marked("", xfine, idxMarked)
        count += 1
    end 

    return xfine, count
end 

calc()
print("calc timing: \n")
@time calc()
xfine, count = calc()
print("length(xfine) = $(length(xfine))\n")
print("count = $count\n")

# d) 
print("d)\n")
x_arr = LinRange(0, 1, 1_000)
y_arr = u.(x_arr)

plot(x_arr, y_arr, title="u", label="u")
savefig("exercise_1_6_d_u.png")

y_arr_initial = u_I.(x_arr, [xcoarse_initial])
y_arr_fine = u_I.(x_arr, [xfine])

plot(x_arr, [y_arr_initial, y_arr_fine], 
    title = "Initial vs. xfine", 
    label=["initial" "fine"]
)
savefig("exercise_1_6_d_initial_vs_xfine.png")

histogram(xfine, title = "Histogram of xfine distribution", label="xfine")
savefig("exercise_1_6_xfine_histogram.png")