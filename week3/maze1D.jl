using Plots
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Random

include("plotpath_2d.jl")
include("plotpath_3d.jl")
include("EToV_50x50_from_image.jl")
include("plotpath.jl")
include("create_maze_EToV.jl")

### Inputs
start, stop, elements = 25, 2476, 2500
cols = 50
c = 1
d = 0

# start, stop, elements = 2, 8, 9 
# EToV = [
#     [1, 2], 
#     [2, 3], [2, 5], 
#     [4, 5], [4, 7], 
#     [5, 6], 
#     [6, 9], 
#     [7, 8], 
#     [8, 9]
# ]

function CO2eq_calculations(CPUtime)
    # PARAMETERS FOR SUUSTAINABILITY CALCULATION
    CO2intensity = 0.100  # [kg CO2/kWh], https://communitiesforfuture.org/collaborate/electricity-map/
    PowerEstimate = 60  # [kW]

    CO2eq = CPUtime / 3600 * PowerEstimate / 1000 * CO2intensity

    return CO2eq
end

function assemble_K(EToV, c, d, start, stop, n_elements)
    E = length(EToV)
    b = zeros(n_elements)

    I = Int[]
    J = Int[]
    V = Float64[]

    for e in 1:E
        i, j = EToV[e]

        push!(I, i); push!(J, i); push!(V,  1)
        push!(I, j); push!(J, j); push!(V,  1)
        push!(I, i); push!(J, j); push!(V, -1)
        push!(I, j); push!(J, i); push!(V, -1)
    end
    A = sparse(I, J, V, n_elements, n_elements)

    # Algorithm 2
    b[start] = c
    b[stop] = d
    
    # Apply boundary conditions
    A[start, :] .= 0
    A[stop, :] .= 0
    A[start, start] = 1
    A[stop, stop] = 1
    
    #don't explicity write zeros
    A = dropzeros(A)    
    return A, b, A \ b 
end
# First compilation 
assemble_K(EToV, c, d,  start, stop, elements)

# A,b = (assemble_K(EToV, c, d,  start, stop, elements))
# @time u = A \ b
# @time A \ b

# @time Array(A) \ b
# plotting_path2d(EToV, u, cols, start)
# plotting_path_3d(EToV, u, cols, start;  zscale = Float64(10))
# plotting_path_3d_4views(EToV, u, cols, start;  zscale = Float64(10))
# plotting_path(EToV, u, cols, start)

CPUtime_list = []
n_list = 3:15
for n ∈ n_list 
    rows, cols = n, n 
    EToV = maze_EToV(rows, cols; seed=42, braid=0.05)   # braid optional
    
    elements = rows * cols
    start = 1
    stop  = rows * cols 
    c, d = 1, 0

    A, b, u = assemble_K(EToV, c, d, start, stop, elements)

    total_time = 0.0 
    fac = 10 
    for i ∈ 1:fac 
        t = @elapsed begin 
            assemble_K(EToV, c, d, start, stop, elements)
        end 
        total_time += t 
    end 
    CPUtime = total_time / fac 
    push!(CPUtime_list, CPUtime)
    
    plotting_path2d(EToV, u, cols, start)
    plotting_path_3d(EToV, u, cols, start; zscale=10.0)
    plotting_path_3d_4views(EToV, u, cols, start; zscale=10.0)
    plotting_path(EToV, u, cols, start)
end 

CO2eq_list = CO2eq_calculations.(CPUtime_list)
println("CPUtime_list = $CPUtime_list")
println("CO2eq = $CO2eq_list")

Plots.plot(n_list, CPUtime_list)
Plots.title!("CPU time vs n × n maze")
Plots.savefig("CPUtime.png")

Plots.plot(n_list, CO2eq_list)
Plots.title!("CO2eq vs n × n maze")
Plots.savefig("CO2eq.png")