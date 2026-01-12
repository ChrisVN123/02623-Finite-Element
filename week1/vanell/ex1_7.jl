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

x_k = LinRange(0, 1, 100)
m_vals = u_I.(x_k,[[0.2,0.3]])
plot(x_k,m_vals)

savefig("cvn1.png")
print(m_vals)

# a, b = 0.0, 1.0
# M = 5
# VX = collect(range(a, b; length=M))          # nodal points (collect is convenient)
# print(VX)
# xs = range(a, b; length=4000)                # fine grid for plotting

# # --- evaluate both ---
# ya = u.(xs)
# yi = [u_I(x, VX) for x in xs]                # u_I is scalar in x, so use comprehension
# print(yi)
# # --- plot overlay ---
# # p = plot(xs, ya, label="Analytical u(x)", linewidth=2, xlabel="x", ylabel="value")
# # plot!(p, xs, yi, label="Interpolant u_I(x)", linewidth=2, linestyle=:dash)

# # # (optional) show nodes and nodal values
# # scatter!(p, VX, u.(VX), label="Nodes", markersize=4)
# # savefig("cvn_test.png")
