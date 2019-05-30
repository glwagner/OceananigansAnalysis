using Oceananigans.Operators

havg(ϕ) = mean(ϕ, dims=(1, 2))
havg(ϕ::Field) = mean(ϕ.data, dims=(1, 2))
means(vars...) = (havg(v.data) for v in vars)

maximum(ϕ::Field) = maximum(ϕ.data)
minimum(ϕ::Field) = minimum(ϕ.data)
abs(ϕ::Field) = abs.(ϕ.data)
absmax(ϕ::Field) = maximum(abs(ϕ))

function fluctuations(vars...)
    types = [typeof(v) for v in vars]
    datas = [v.data .- havg(v) for v in vars]
    return (T(datas[i], vars[i].grid) for (i, T) in enumerate(types))
end

function fluctuation(v)
    Field = typeof(v)
    data = v.data .- havg(v)
    return Field(data, v.grid)
end

function normalize!(ϕ::Field)
    ϕ.data .-= minimum(mean(ϕ.data, dims=(1, 2)))
    ϕ.data ./= maximum(mean(ϕ.data, dims=(1, 2)))
    return nothing
end

@inline incmod1(a, n) = ifelse(a==n, 1, a + 1)
@inline decmod1(a, n) = ifelse(a==1, n, a - 1)

#
# Derivatives... easy case: parallel strain
#

∂x(u::FaceFieldX, i, j, k) = δx_f2c(u.grid, u, i, j, k) / u.grid.Δx
∂y(v::FaceFieldY, i, j, k) = δy_f2c(v.grid, v, i, j, k) / v.grid.Δy
∂z(w::FaceFieldZ, i, j, k) = δz_f2c(w.grid, w, i, j, k) / w.grid.Δz

#
# Derivatives: "edge cases" 😏: Vorticity
#

#=
z
^

(i, k-1)   * ----- *  (i+1, k-1)
           |       |
           u       |
           |       |
 (i, k)    * --w-- *   (i+1, k)

=#

function ∂x(w::FaceFieldZ, i, j, k)
    if k == 1
        return 0.25/w.grid.Δx * (δx_f2e(w.grid, w, i, j, k)   + δx_f2e(w.grid, w, incmod1(i, w.grid.Nz), j, k))
    else
        return 0.25/w.grid.Δx * (   δx_f2e(w.grid, w, i, j, k)   + δx_f2e(w.grid, w, incmod1(i, w.grid.Nz), j, k)
                                  + δx_f2e(w.grid, w, i, j, k-1) + δx_f2e(w.grid, w, incmod1(i, w.grid.Nz), j, k-1) )
    end
end

function ∂z(u::FaceFieldX, i, j, k)
    if k == 1
        return 0.25/u.grid.Δz * (δz_f2e(u.grid, u, i, j, k)   + δz_f2e(u.grid, u, incmod1(i, u.grid.Nx), j, k))
    else
        return 0.25/u.grid.Δz * (   δz_f2e(u.grid, u, i, j, k)   + δz_f2e(u.grid, u, incmod1(i, u.grid.Nx), j, k)
                                  + δz_f2e(u.grid, u, i, j, k-1) + δz_f2e(u.grid, u, incmod1(i, u.grid.Nx), j, k-1) )
    end
end

#=
z
^

(j, k-1)   * ----- *  (j+1, k-1)
           |       |
           v       |
           |       |
 (j, k)    * --w-- *   (j+1, k)

=#

function ∂y(w::FaceFieldZ, i, j, k)
    if k == 1
        return 0.25/w.grid.Δy * (δy_f2e(w.grid, w, i, j, k)   + δy_f2e(w.grid, w, i, incmod1(j, w.grid.Ny), k))
    else
        return 0.25/w.grid.Δy * (   δy_f2e(w.grid, w, i, j, k)   + δy_f2e(w.grid, w, i, incmod1(j, w.grid.Ny), k)
                                  + δy_f2e(w.grid, w, i, j, k-1) + δy_f2e(w.grid, w, i, incmod1(j, w.grid.Ny), k-1) )
    end
end

function ∂z(v::FaceFieldY, i, j, k)
    if k == 1 # surface
        return 0.25/v.grid.Δz * (δz_f2e(v.grid, v, i, j, k) + δz_f2e(v.grid, v, i, incmod1(j, v.grid.Ny), k))
    else
        return 0.25/v.grid.Δz * (   δz_f2e(v.grid, v, i, j, k)   + δz_f2e(v.grid, v, i, incmod1(j, v.grid.Ny), k)
                                  + δz_f2e(v.grid, v, i, j, k-1) + δz_f2e(v.grid, v, i, incmod1(j, v.grid.Ny), k-1) )
    end
end


#=
y
^

(i, j+1)   * ----- *  (i+1, j+1)
           |       |
           u       |
           |       |
 (i, j)    * --v-- *   (i+1, j)

=#

function ∂y(u::FaceFieldX, i, j, k)
    return 0.25/u.grid.Δy * (   δy_f2e(u.grid, u, i, j, k)   + δy_f2e(u.grid, u, incmod1(i, u.grid.Nx), j, k)
                              + δy_f2e(u.grid, u, i, incmod1(j, u.grid.Ny), k) + δy_f2e(u.grid, u, incmod1(i, u.grid.Nx), incmod1(j, u.grid.Ny), k) )
end

function ∂x(v::FaceFieldY, i, j, k)
    return 0.25/v.grid.Δx * (   δx_f2e(v.grid, v, i, j, k)   + δx_f2e(v.grid, v, incmod1(i, v.grid.Nx), j, k)
                              + δx_f2e(v.grid, v, i, incmod1(j, v.grid.Ny), k) + δx_f2e(v.grid, v, incmod1(i, v.grid.Nx), incmod1(j, v.grid.Ny), k) )
end

function kinetic_energy(u, v, w)
    edata = zeros(size(u))
    nx, ny, nz = size(u)

    for k = 1:nz, j = 1:ny, i = 1:nx
        @inbounds begin
            edata[i, j, k] = 0.25 * (  (u[i, j, k] + u[incmod1(i, nx), j, k])^2
                                     + (v[i, j, k] + v[i, incmod1(j, ny), k])^2)
            if k < nz
                edata[i, j, k] += 0.25 * (w[i, j, k] + w[i, j, k+1])^2
            else
                edata[i, j, k] += 0.25 * w[i, j, k]^2
            end
        end
    end

    return CellField(edata, u.grid)
end

function dissipation(ν, u, v, w)
    ϵ = zeros(size(u))
    nx, ny, nz = size(u)

    for k = 1:nz, j = 1:ny, i = 1:nx
        @inbounds begin
            ϵ[i, j, k] = ν * (
                  ∂x(u, i, j, k)^2 + ∂y(u, i, j, k)^2 + ∂z(u, i, j, k)^2
                + ∂x(v, i, j, k)^2 + ∂y(v, i, j, k)^2 + ∂z(v, i, j, k)^2
                + ∂x(w, i, j, k)^2 + ∂y(w, i, j, k)^2 + ∂z(w, i, j, k)^2 )
        end
    end

    return CellField(ϵ, u.grid)
end

function *(w::FaceFieldZ, c::CellField)
    data = zeros(size(w))
    nx, ny, nz = size(w)

    for k = 1:nz, j = 1:ny, i = 1:nx
        if k < nz
            @inbounds data[i, j, k] = 0.5*(w[i, j, k] + w[i, j, k+1]) * c[i, j, k]
        else # Use no-penetration condition which implies w=0 at bottom.
            @inbounds data[i, j, k] = 0.5*w[i, j, k]*c[i, j, k]
        end
    end

    return CellField(data, c.grid)
end
