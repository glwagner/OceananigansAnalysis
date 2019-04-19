import PyPlot: plot

xnodes(ϕ::Field) = reshape(ϕ.grid.xC, ϕ.grid.Nx, 1, 1)
ynodes(ϕ::Field) = reshape(ϕ.grid.yC, 1, ϕ.grid.Ny, 1)
znodes(ϕ::Field) = reshape(ϕ.grid.zC, 1, 1, ϕ.grid.Nz)

xnodes(ϕ::FaceFieldX) = reshape(ϕ.grid.xF[1:end-1], ϕ.grid.Nx, 1, 1)
ynodes(ϕ::FaceFieldY) = reshape(ϕ.grid.yF[1:end-1], 1, ϕ.grid.Ny, 1)
znodes(ϕ::FaceFieldZ) = reshape(ϕ.grid.zF[1:end-1], 1, 1, ϕ.grid.Nz)

xzsliceplot(ϕ::Field, slice=1; kwargs...) = pcolormesh(
    repeat(dropdims(xnodes(ϕ), dims=2), 1, ϕ.grid.Nz),
    repeat(dropdims(znodes(ϕ), dims=2), ϕ.grid.Nx, 1),
    ϕ.data[slice, :, :]; kwargs...)

aspectratio(a, ax=gca(); adjustable="box") = ax.set_aspect(a, adjustable=adjustable)
makesquare(ax=gca()) = aspectratio(1, ax)
makesquare(axs::AbstractArray) = for ax in axs; makesquare(ax); end

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
