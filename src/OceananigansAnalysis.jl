module OceananigansAnalysis

export parse_filename, num, name_without_num, load_grid, load_solution, 
       load_solution!, sort_paths, simtime, simiter, augment_vars!,

       havg, means, fluctuations, maximum, minimum,

       makesquare, xzsliceplot, usecmbright, removespines, cornerspines,

       create_timeseries, iterations, times, getdata, getkappa, get_data_params

using NetCDF, Glob, PyPlot, Oceananigans, Statistics, JLD2

import Base: maximum, minimum

# ρ₀ = 1027 kg/m³, cp = 4181.3 J/(kg·K), α = 2.07e-4 K⁻¹, g = 9.80665 m/s²
const ρ₀ = 1027.0
const cP = 4181.3
const α = 2.07e-4
const g = 9.80665
const f = 1e-4

include("plotting.jl")
include("parsing.jl")
include("timeseries_analysis.jl")

buoyancy_flux(Q) = -α * g * Q / (ρ₀ * cP)
velocity_flux(wind_stress) = wind_stress / ρ₀
buoyancy_gradient(dTdz) = α * g * dTdz

havg(ϕ) = mean(ϕ, dims=(1, 2))
havg(ϕ::Field) = mean(ϕ.data, dims=(1, 2))
means(vars...) = (havg(v.data) for v in vars)

maximum(ϕ::Field) = maximum(ϕ.data)
minimum(ϕ::Field) = minimum(ϕ.data)

function fluctuations(vars...)
    types = [typeof(v) for v in vars]
    datas = [v.data .- havg(v) for v in vars]
    return (T(datas[i], vars[i].grid) for (i, T) in enumerate(types))
end

#=
function turbulent_kinetic_energy(u, v, w)
    nx, ny, nz = size(u)
    e = zeros(nx, ny, nz)

    for k in 1:nz, j=1:ny, i=1:nx
        e[i, j, k] = avgxc2f(u.grid, u.data, i, j, k)^2 

end
=#

function create_timeseries(timeseriespath; simname="", dir=".", noutput=nothing, verbose=false)
    datapaths = sort_paths(glob(simname * "*.nc", dir))

    if verbose
        @info "Creating a timeseries from the following data:\n $(("$d \n" for d in datapaths)...)"
    end

    if noutput == nothing
        noutput = length(datapaths) # not necessarily reliable...
    end

    # Save simulation metadata
    sample_datapath = datapaths[1]
    grid = load_grid(sample_datapath)
    vars = parse_filename(sample_datapath, ["N", "tau", "Q", "dTdz", "k", "dt", "days"])
    augment_vars!(vars)

    jldopen(timeseriespath, "a+") do file
        file["grid/N"] = grid.Nz
        file["grid/L"] = grid.Lz
    end

    jldopen(timeseriespath, "a+") do file
        file["boundary_conditions/Fb"] = vars["Fb"] 
        file["boundary_conditions/Fu"] = vars["Fu"] 
        file["initial_condition/Bz"] = vars["dbdz"]

        file["timestepping/dt"] = vars["dt"]
        file["timestepping/tfinal"] = vars["tfinal"]

        file["constants/ρ₀"] = ρ₀
        file["constants/cP"] = cP
        file["constants/g"] = g
        file["constants/α"] = α
        file["constants/f"] = f
        file["constants/κ"] = vars["k"]
    end

    # Initialize.
    u, v, w, θ, s = load_solution(datapaths[1])

    for (idata, datapath) in enumerate(datapaths)
        iter = simiter(datapath, noutput; prefix=simname)
           t = simtime(datapath, noutput; prefix=simname)

        if verbose
            @info "Processing iteration $iter"
            t0 = time_ns()
        end

        # Conserve memory, or attempt to.
        load_solution!(u, v, w, θ, s, datapath)

        U, V, W, T, S = means(u, v, w, θ, s)

        if verbose
            @info "... processing took $((time_ns()-t0)) s. Saving..."
            t0 = time_ns()
        end

        jldopen(timeseriespath, "a+") do file
            file["timeseries/t/$iter"] = t
            file["timeseries/U/$iter"] = U
            file["timeseries/V/$iter"] = V
            file["timeseries/T/$iter"] = T
            file["timeseries/S/$iter"] = S
        end

        if verbose
            @info "... and saving took $((time_ns()-t0)) s."
        end
    end

    return nothing
end

end # module
