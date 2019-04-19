module OceananigansAnalysis

export  parse_filename, load_grid, load_solution, simtime, simiter,
        havg, means, fluctuations, maximum, minimum,
        makesquare, xzsliceplot,
        create_timeseries

using NetCDF, Glob, PyPlot, Oceananigans, Statistics, JLD2

import Base: maximum, minimum

# ρ₀ = 1027 kg/m³, cp = 4181.3 J/(kg·K), α = 2.07e-4 K⁻¹, g = 9.80665 m/s²
const ρ₀ = 1027.0
const cP = 4181.3
const α = 2.07e-4
const g = 9.80665

include("plotting.jl")
include("parsing.jl")

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

function create_timeseries(timeseriespath; simname="", dir=".", noutput=nothing)
    datapaths = glob(simname * "*.nc", dir)

    if noutput == nothing
        noutput = length(datapaths) # not necessarily reliable...
    end

    grid = load_grid(datapaths[1])

    jldopen(timeseriespath, "a+") do file
        file["grid/N"] = grid.Nz
        file["grid/L"] = grid.Lz
    end

    for datapath in datapaths
        iter = simiter(simname, datapath, noutput)

        u, v, w, θ, s = load_solution(datapath)
        U, V, W, T, S = means(u, v, w, θ, s)

        jldopen(timeseriespath, "a+") do file
            file["timeseries/t/$iter"] = simtime(simname, datapath, noutput)
            file["timeseries/U/$iter"] = U
            file["timeseries/V/$iter"] = V
            file["timeseries/T/$iter"] = T
            file["timeseries/S/$iter"] = S
        end
    end

    return nothing
end

end # module
