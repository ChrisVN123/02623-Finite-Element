using Plots 
using BenchmarkTools

function f(x)
    return 2*x^2+x-1
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

function u_I(x, node, VX)
    res = 0
    for i in x
        res += N(x, node, VX)
    end
    return res
end

VX = [1,2,3,4]
print(Ref(VX))
#print(N(1,1,VX))
x_vals = 1:0.1:4
fun = u_I.(x_vals, 2, Ref(VX))
fun2 = u_I.(x_vals, 1, Ref(VX))
plot(x_vals, fun)
plot!(x_vals,fun2)
savefig("test.png")


# test = N(0.5, 1, 3, 3)
# print(test)


# test = [1, 2, 3, 4]
# let resu = 0
#     for x in test
#         resu += x
#     end
# end
# print(resu)


