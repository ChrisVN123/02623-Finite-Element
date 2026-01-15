using Printf
using PyCall
using PyPlot

include("driver28b.jl")
include("driver28c.jl")

# TODO: PUT IN YOUR GROUP NO AND STUDENT IDs FOR THE GROUP HERE
groupNo = "X"  # Replace "X" with your group no.
groupStudentIDs = "Y/Z"  # Replace "Y/Z" with your student id's.

# PATHS
dirThisScript = @__DIR__
dirStoreResults = "/Users/apek/02623/Handinresults/"

# PARAMETERS FOR SUSTAINABILITY CALCULATION
CO2intensity = 0.285  # [kg CO2/kWh]
PowerEstimate = 60  # [kW]

# PARAMETERS FOR THE SELECTED EXERCISES (DO NOT CHANGE)
x0, y0 = 0, 0
L1, L2 = 1, 1
noelms1, noelms2 = 40, 50
lam1, lam2 = 1, 1
fun = (x, y) -> cos(pi * x) * cos(pi * y)
qt = (x, y) -> 2 * pi^2 * cos(pi * x) * cos(pi * y)

# EXECUTE CODE
# Call Group 30 solver
start_time = time()
VX, VY, EToV, U = Driver28b(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, fun, qt)
tend = time() - start_time
DOF1 = length(vec(U))

x0, y0 = -1, -1
L1, L2 = 2, 2
start_time = time()
VX2, VY2, EToV2, U2 = Driver28c(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, fun, qt)
tend2 = time() - start_time
DOF2 = length(vec(U2))

CPUtime1 = tend
CPUtime2 = tend2
CO2eq1 = CPUtime1 / 3600 * PowerEstimate / 1000 * CO2intensity
CO2eq2 = CPUtime2 / 3600 * PowerEstimate / 1000 * CO2intensity

# Visualization
tri = pyimport("matplotlib.tri")
fig, ax = subplots(1, 2, figsize=(15, 6))
# Matplotlib expects 0-based triangle indices.
triang = tri.Triangulation(VX, VY, EToV .- 1)
ax[1].tripcolor(triang, U, shading="flat")
ax[1].set_title(@sprintf(
    "2.8b. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e",
    groupStudentIDs,
    tend,
    DOF1,
    noelms1,
    noelms2,
    CO2eq1
))

triang2 = tri.Triangulation(VX2, VY2, EToV2 .- 1)
ax[2].tripcolor(triang2, U2, shading="flat")
ax[2].set_title(@sprintf(
    "2.8c. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e",
    groupStudentIDs,
    tend2,
    DOF2,
    noelms1,
    noelms2,
    CO2eq2
))

show()

# STORE THE RESULTS
filename = @sprintf("Week2ResultsGroup%s.txt", groupNo)
open(joinpath(dirStoreResults, filename), "w") do file
    write(file, @sprintf("%s\n", groupNo))
    write(file, @sprintf("%s\n", groupStudentIDs))
    write(file, @sprintf("%.4e %d %.4e\n", CPUtime1, DOF1, CO2eq1))
    write(file, @sprintf("%.4e %d %.4e\n", CPUtime2, DOF2, CO2eq2))
end
