module OceananigansAnalysis

export # misc constants
    α, g,

    # file wrangling
    parse_filename, num, name_without_num, load_grid, load_solution,
    load_solution!, sort_paths, simtime, simiter, augment_vars!,
    diffusivity, create_timeseries, iterations, times, getdata, getkappa,
    get_data_params, makesquare, makesimple,

    # operations
    logabs, maxabs, havg, means, fluctuation, fluctuations, normalize!,
    kinetic_energy, turbulent_kinetic_energy, dissipation, NewField,
    Umax, Δmin, *,

    # plotting
    makesquare, plot_xzslice, plot_hmean, usecmbright, removespines, cornerspines, plot_havg,
    three_field_viz, allsides

using NetCDF, Glob, PyPlot, Oceananigans, Statistics, JLD2

import Base: maximum, minimum, abs, *, /, log

zerofunk(args...) = 0

function set_ic!(model; ics...)
    for (fld, ic) in ics
        if fld ∈ (:u, :v, :w)
            ϕ = getproperty(model.velocities, fld)
        else
            ϕ = getproperty(model.tracers, fld)
        end
        data(ϕ) .= ic.(nodes(ϕ)...)
    end
    return nothing
end

# ρ₀ = 1027 kg/m³, cp = 4181.3 J/(kg·K), α = 2.07e-4 K⁻¹, g = 9.80665 m/s²
const ρ₀ = 1027.0
const cP = 4181.3
const α = 2.07e-4
const g = 9.80665
const f = 1e-4

function dotnest(f, a)
    undotted = Expr(:call, f, a)
    return :(@. $undotted)
end

function NewField(c::Field, fns...)
    ex = Expr(:call, fns[1], c.data)
    ex = :(@. $ex)
    for fn in fns
        ex = Expr(:call, fn, ex)
        ex = :(@. $ex)
    end
    eval(:(data = @. $ex))
    return typeof(c)(data, c.grid)
end

log(c::CellField) = CellField(log.(c.data), c.grid)
logabs(c::CellField) = CellField(log.(abs.(c.data)), c.grid)

include("plotting.jl")
include("file_wrangling.jl")
include("timeseries_analysis.jl")
include("operators.jl")

buoyancy_flux(Q) = -α * g * Q / (ρ₀ * cP)
velocity_flux(wind_stress) = wind_stress / ρ₀
buoyancy_gradient(dTdz) = α * g * dTdz

end # module
