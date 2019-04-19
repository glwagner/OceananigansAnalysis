function parse_filename(simname, filepath)
    filename = basename(filepath)
    filename = filename[length(simname)+2:end]

    varnames = ["N", "Q", "dTdz", "k", "dt", "days"]

    if simname == "wind_stress"
        push!("tau", varnames)
    end

    vartypes = cat([Int], [Float64 for i = 1:length(varnames)-1], dims=1)
    vars = Dict{String, Any}()

    for (i, v) in enumerate(varnames)
        start = findnext(v, filename, 1)[end] + 1
        finish = findnext("_", filename, start)[end] - 1
        vars[v] = parse.(vartypes[i], filename[start:finish])
    end

    daysstart = findnext("days", filename, 1)[end] + 1
    numstart = findnext("_", filename, daysstart)[end] + 1
    numfinish = findnext(".", filename, numstart)[end] - 1

    vars["num"] = parse(Int, filename[numstart:numfinish])
    vars["Fb"] = buoyancy_flux(vars["Q"])
    vars["dbdz"] = buoyancy_gradient(vars["dTdz"])
    vars["tfinal"] = vars["days"] * 24 * 3600

    return vars
end

function simtime(simname, filepath, noutput)
    vars = parse_filename(simname, filepath)
    return vars["num"]/(noutput-1) * vars["tfinal"]
end

function iter(simname, filepath, noutput)
    vars = parse_filename(simname, filepath)
    return Int(vars["num"]/(noutput-1) * vars["tfinal"] / vars["dt"])
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

