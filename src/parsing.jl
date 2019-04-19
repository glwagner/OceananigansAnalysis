function parse_filename(filepath, varnames, vartypes=[Float64 for i=1:length(varnames)]; prefix=nothing)
    filename = basename(filepath)
    
    if prefix == nothing
        iprefix = 1
    else
        prefixrange = findnext(prefix, filename, 1)
        prefixrange == nothing && error("prefix $prefix not found in 
                                        basename $filename associated with path $filepath")
        iprefix = prefixrange[end]
    end

    vars = Dict{String, Any}()

    for (i, v) in enumerate(varnames)
        start = findnext(v, filename, iprefix)[end] + 1
        finish = findnext("_", filename, start)[end] - 1
        vars[v] = parse.(vartypes[i], filename[start:finish])
    end

    vars["num"] = num(filepath)

    return vars
end

"""
    num(filepath)

Return a number appended to the end of a filepath, such that

<stuff>_num.<file_suffix>
"""
function num(filepath)
    filename = basename(filepath)
    numrange = findprev("_", filename, length(filename))

    numrange == nothing && error("No underscore was found in $filepath")

    numstart = numrange[end] + 1
    numfinish = findnext(".", filename, numstart)[end] - 1
    return parse(Int, filename[numstart:numfinish])
end

function name_without_num(filepath)
    filename = basename(filepath)
    numrange = findprev("_", filename, length(filename))
    beforenum = numrange[1] - 1
    numrange == nothing && error("No underscore was found in $filepath")
    return filename[1:beforenum]
end

function sort_paths(filepaths)
    nums = num.(filepaths)
    ii = sortperm(nums)
    return filepaths[ii]
end

function augment_vars!(vars)
    vars["N"] = Int(vars["N"])
    vars["Fu"] = velocity_flux(vars["tau"])
    vars["Fb"] = buoyancy_flux(vars["Q"])
    vars["dbdz"] = buoyancy_gradient(vars["dTdz"])
    vars["tfinal"] = vars["days"] * 24 * 3600
    return nothing
end

function simtime(filepath, noutput::Int; prefix=nothing)
    vars = parse_filename(filepath, ["days"]; prefix=prefix)
    tfinal = vars["days"] * 24 * 3600
    return vars["num"]/noutput * tfinal
end

function simiter(filepath, noutput::Int; prefix=nothing)
    vars = parse_filename(filepath, ["days", "dt"]; prefix=prefix)
    tfinal = vars["days"] * 24 * 3600
    return Int(vars["num"]/noutput * tfinal / vars["dt"])
end


function load_grid(filepath)
    xC = ncread(filepath, "xC")
    yC = ncread(filepath, "yC")
    zC = ncread(filepath, "zC")

    xF = ncread(filepath, "xF")
    yF = ncread(filepath, "yF")
    zF = ncread(filepath, "zF")

    Nx = length(xC)
    Ny = length(yC)
    Nz = length(zC)

    Δx = xC[2] - xC[1]
    Δy = yC[2] - yC[1]
    Δz = zC[1] - zC[2]

    Lx = Δx * Nx
    Ly = Δy * Ny
    Lz = Δz * Nz

    return RegularCartesianGrid((Nx, Ny, Nz), (Lx, Ly, Lz))
end


function load_solution(filepath)
    grid = load_grid(filepath)

    u = FaceFieldX(ncread(filepath, "u"), grid)
    v = FaceFieldY(ncread(filepath, "v"), grid)
    w = FaceFieldZ(ncread(filepath, "w"), grid)

    θ = CellField(ncread(filepath, "T"), grid)
    s = CellField(ncread(filepath, "S"), grid)

    return u, v, w, θ, s
end

