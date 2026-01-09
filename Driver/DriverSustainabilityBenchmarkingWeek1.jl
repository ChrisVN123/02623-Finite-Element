# using Printf
# using PyPlot
using Plots

# Import the DriverAMR17 function from the external Julia script DriverAMR17.jl
include("DriverAMR17.jl")

# Define input parameters (DO NOT CHANGE THIS PART)
funu(x) = exp(-800 * (x - 0.4)^2) + 0.25 * exp(-40 * (x - 0.8)^2)
func(x) = (-1601 * exp(-800 * (x - 0.4)^2) + (-1600 * x + 640.0)^2 * exp(-800 * (x - 0.4)^2) - 20.25 * exp(-40 * (x - 0.8)^2) + 0.25 * (-80 * x + 64.0)^2 * exp(-40 * (x - 0.8)^2))

x = [0.0, 0.5, 1.0] # Initial mesh configuration - do not change
M = length(x)
L = 1
c = funu(x[1])
d = funu(x[end])
tol = 1e-4
maxit = 50

# Let's call the FEM BVP 1D Solver with AMR
# time the code using time
fac = 100 # we do multiple runs to get the average time
start_time = time()

xAMR = 0
u = 0
iter = 0

for i in 1:fac
    global xAMR, u, iter = DriverAMR17(L, c, d, x, func, tol, maxit)
end
CPUtime = (time() - start_time) / fac


# Plot
DOF = length(xAMR)
CO2eq = CPUtime / 3600 * 86 / 1000 * 0.135 # valid for MacBook Pro (assumed power consumption 105)
print("CO2eq = $CO2eq\n")
plot(
    xAMR, funu.(xAMR), 
    title="Group: <Your Group ID>, Iter: $(iter), Time: $(round(CPUtime, digits=4)) s, DOF: $(DOF), CO2e=$(CO2eq) kg CO2",
    label=["Exact", "AMR"]
)
plot!(xAMR, u)
savefig("DriverSustainability.png")

# Element size distribution
h = diff(xAMR)
# figure()
histogram(h)
savefig("DriverSustainabilityHistogram.png")
# xlabel("h")
# ylabel("# of elements")
# show()
