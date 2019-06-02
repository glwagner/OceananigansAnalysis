import PyPlot: plot

const allsides = ("top", "bottom", "left", "right")

#
# Helper functions
#

xnodes(ϕ::Field) = reshape(ϕ.grid.xC, ϕ.grid.Nx, 1, 1)
ynodes(ϕ::Field) = reshape(ϕ.grid.yC, 1, ϕ.grid.Ny, 1)
znodes(ϕ::Field) = reshape(ϕ.grid.zC, 1, 1, ϕ.grid.Nz)

xnodes(ϕ::FaceFieldX) = reshape(ϕ.grid.xF[1:end-1], ϕ.grid.Nx, 1, 1)
ynodes(ϕ::FaceFieldY) = reshape(ϕ.grid.yF[1:end-1], 1, ϕ.grid.Ny, 1)
znodes(ϕ::FaceFieldZ) = reshape(ϕ.grid.zF[1:end-1], 1, 1, ϕ.grid.Nz)

nodes(ϕ) = (xnodes(ϕ), ynodes(ϕ), znodes(ϕ))

zslice(ϕ) = repeat(dropdims(znodes(ϕ), dims=2), ϕ.grid.Nx, 1)
xslice(ϕ) = repeat(dropdims(xnodes(ϕ), dims=2), 1, ϕ.grid.Nz)

x3d(ϕ) = repeat(xnodes(ϕ), 1, ϕ.grid.Ny, ϕ.grid.Nz)
y3d(ϕ) = repeat(ynodes(ϕ), ϕ.grid.Nx, 1, ϕ.grid.Nz)
z3d(ϕ) = repeat(znodes(ϕ), ϕ.grid.Nx, ϕ.grid.Ny, 1)

#
# Plot types
#

plot_xzslice(ϕ, slice=1, args...; kwargs...) = 
    pcolormesh(view(x3d(ϕ), :, slice, :), view(z3d(ϕ), :, slice, :), 
               view(data(ϕ), :, slice, :), args...; kwargs...)

plot_xyslice(ϕ, slice=1, args...; kwargs...) = 
    pcolormesh(view(x3d(ϕ), :, :, slice), view(y3d(ϕ), :, :, slice), 
               view(data(ϕ), :, :, slice), args...; kwargs...)

function plot_hmean(f::Function, ϕ, args...; normalize=false, kwargs...)
    ϕhmean = dropdims(mean(data(ϕ), dims=(1, 2)), dims=(1, 2))
    if !normalize
        ϕnorm = 1
    else
        ϕmean = mean(data(ϕ))
        ϕhmean = ϕhmean .- ϕmean
        ϕnorm = maximum(ϕhmean) - minimum(ϕhmean)
    end
    PyPlot.plot(f.(ϕhmean/ϕnorm), ϕ.grid.zC, args...; kwargs...)
end

plot_hmean(ϕ::Field, args...; kwargs...) =
    plot_hmean(x->x, ϕ, args...; kwargs...)

#
# Plot modifiers
#

aspectratio(a, ax=gca(); adjustable="box") = ax.set_aspect(a, adjustable=adjustable)
makesquare(ax=gca()) = aspectratio(1, ax)
makesquare(axs::AbstractArray) = for ax in axs; makesquare(ax); end

function makesimple(axs)
    for ax in axs
        ax.set_aspect(1)
        ax.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false)
        removespines(allsides..., ax=ax)
    end
    return nothing
end

latexpreamble = """
\\usepackage{cmbright}
\\renewcommand{\\b}[1]    {\\boldsymbol{#1}}
\\renewcommand{\\r}[1]    {\\mathrm{#1}}
\\renewcommand{\\d}       {\\partial}
"""

function usecmbright()
  rc("text.latex", preamble=latexpreamble)
  rc("font", family="sans-serif")
  nothing
end

"Remove `spine` from `ax`."
function removespine(side, ax=gca())
    ax.spines[side].set_visible(false)
    keywords = Dict(Symbol(side)=>false, Symbol(:label, side)=>false)
    ax.tick_params(keywords)
    nothing
end

removespines(sides...; ax=gca()) = for side in sides; removespine(side, ax); end
cornerspines() = removespines("top", "right")

function three_field_viz(fields...; seq=(false, true), clips=(0.5, 0.5))

    fig, axs = subplots(ncols=3, figsize=(12, 4))

    allsides = ("top", "bottom", "left", "right")

    for ax in axs
        ax.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false)
        sca(ax)
        removespines(allsides...)
    end

    [makesquare(axs[i]) for i=1:2]

    lims = []
    cmaps = []
    for (i, s) in enumerate(seq)
        if s
            push!(lims, [0, 1]*clips[i]*absmax(fields[i]))
            push!(cmaps, "YlGnBu_r")

        else
            push!(lims, [-1, 1]*clips[i]*absmax(fields[i]))
            push!(cmaps, "RdBu_r")
        end
    end

    axs[1].tick_params(left=true, labelleft=true)
    axs[3].tick_params(right=true, labelright=true)

    for (i, ϕ) in enumerate(fields)
        sca(axs[i])
        xzsliceplot(ϕ, cmap=cmaps[i], vmin=lims[i][1], vmax=lims[i][2])
        normalize!(ϕ)
    end

    axs[3].spines["bottom"].set_visible(true)
    axs[3].spines["right"].set_visible(true)
    axs[3].tick_params(labelbottom=true)

    sca(axs[3])
    for ϕ in fields
        plot_havg(ϕ)
    end

    return fig, axs
end

function fluctuations_and_means(field1, field2, meanfields...; seq=(false, true),
                                clipping=(0.5, 0.5), mean_names=["" for ϕ in meanfields])

    fig, axs = subplots(ncols=3, figsize=(12, 4))

    allsides = ("top", "bottom", "left", "right")

    for ax in axs
        ax.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false)
        sca(ax)
        removespines(allsides...)
    end

    [makesquare(axs[i]) for i=1:2]

    fields = [field1, field2]
    lims, cmaps = [], []
    for (i, s) in enumerate(seq)
        if s
            push!(lims, [0, 1]*clips[i]*absmax(fields[i]))
            push!(cmaps, "YlGnBu_r")
        else
            push!(lims, [-1, 1]*clips[i]*absmax(fields[i]))
            push!(cmaps, "RdBu_r")
        end
    end

    axs[1].tick_params(left=true, labelleft=true)
    axs[3].tick_params(right=true, labelright=true)

    for (i, ϕ) in enumerate(fields)
        sca(axs[i])
        xzsliceplot(ϕ, cmap=cmaps[i], vmin=lims[i][1], vmax=lims[i][2])
    end

    axs[3].spines["bottom"].set_visible(true)
    axs[3].spines["right"].set_visible(true)
    axs[3].tick_params(labelbottom=true)

    sca(axs[3])
    for (i, ϕ) in enumerate(meanfields)
        normalize!(ϕ)
        plot_havg(ϕ, label=mean_names[i])
    end

    return fig, axs
end
