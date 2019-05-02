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

function diffusivity(filepath; prefix=nothing)
    vars = parse_filename(filepath, ["k"]; prefix=prefix)
    return vars["k"]
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


function load_solution!(u, v, w, θ, s, filepath)
    grid = load_grid(filepath)

    ncread!(filepath, "u", u.data)
    ncread!(filepath, "v", v.data)
    ncread!(filepath, "w", w.data)
    ncread!(filepath, "T", θ.data)
    ncread!(filepath, "S", s.data)

    return nothing
end

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

