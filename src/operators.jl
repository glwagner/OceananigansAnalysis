using
    Oceananigans,
    Oceananigans.Operators,
    Oceananigans.TurbulenceClosures

using Oceananigans.TurbulenceClosures: ▶x_caa, ▶y_aca, ▶z_aac, ▶x_faa, ▶y_afa, ▶z_aaf

havg(ϕ) = mean(ϕ, dims=(1, 2))
havg(ϕ::Field) = mean(data(ϕ), dims=(1, 2))
means(vars...) = (havg(data(v)) for v in vars)

maximum(ϕ::Field) = maximum(data(ϕ))
minimum(ϕ::Field) = minimum(data(ϕ))
abs(ϕ::Field) = abs.(data(ϕ))
maxabs(ϕ::Field) = maximum(abs(ϕ))

Umax(u, v, w) = max(maxabs(u), maxabs(v), maxabs(w))
Umax(model) = Umax(model.velocities.u, model.velocities.v, model.velocities.w)

Δmin(grid) = min(grid.Δx, grid.Δy, grid.Δz)

fieldtype(::CellField) = CellField
fieldtype(::FaceFieldX) = FaceFieldX
fieldtype(::FaceFieldY) = FaceFieldY
fieldtype(::FaceFieldZ) = FaceFieldZ

function fluctuations(vars...)
    types = [fieldtype(v) for v in vars]
    datas = [v.data .- havg(v) for v in vars]
    return (T(datas[i], vars[i].grid) for (i, T) in enumerate(types))
end

function fluctuation(v)
    Field = fieldtype(v)
    data = v.data .- havg(v)
    return Field(data, v.grid)
end

function normalize!(ϕ::Field)
    ϕ.data .-= minimum(mean(ϕ.data, dims=(1, 2)))
    ϕ.data ./= maximum(mean(ϕ.data, dims=(1, 2)))
    return nothing
end

function *(u::FaceFieldX, ϕ::CellField)
    uϕdata = zeros(size(data(u))...)

    for k = 1:u.grid.Nz, j = 1:u.grid.Ny, i = 1:u.grid.Nx
        @inbounds uϕdata[i, j, k] = ▶x_caa(i, j, k, u.grid, u.data) * ϕ[i, j, k]
    end

    return CellField(uϕdata, ϕ.grid)
end

function *(v::FaceFieldY, ϕ::CellField)
    vϕdata = zeros(size(v)...)

    for k = 1:v.grid.Nz, j = 1:v.grid.Ny, i = 1:v.grid.Nx
        @inbounds vϕdata[i, j, k] = ▶y_aca(i, j, k, v.grid, v.data) * ϕ[i, j, k]
    end

    return CellField(vϕdata, ϕ.grid)
end

function *(w::FaceFieldZ, ϕ::CellField)
    wϕ = CellField(CPU(), ϕ.grid)

    for k = 1:w.grid.Nz, j = 1:w.grid.Ny, i = 1:w.grid.Nx
        @inbounds wϕ[i, j, k] = ▶z_aac(i, j, k, w.grid, w.data) * ϕ[i, j, k]
    end

    return wϕ
end

Base.@propagate_inbounds pow2(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

function kinetic_energy(u, v, w)

    edata = zeros(CPU(), u.grid)

    for k = 1:u.grid.Nz, j = 1:u.grid.Ny, i = 1:u.grid.Nx
        @inbounds edata[i, j, k] = 0.5 * (
              ▶x_caa(i, j, k, u.grid, pow2, u.data)
            + ▶y_aca(i, j, k, u.grid, pow2, v.data)
            + ▶z_aac(i, j, k, u.grid, pow2, w.data)
           )
    end

    return CellField(edata, u.grid)
end

function turbulent_kinetic_energy(u, v, w)
    u′, v′, w′ = fluctuations(u, v, w)
    return kinetic_energy(u′, v′, w′)
end

turbulent_kinetic_energy(model) = turbulent_kinetic_energy(model.velocities.u,
                                                           model.velocities.v,
                                                           model.velocities.w
                                                          )
