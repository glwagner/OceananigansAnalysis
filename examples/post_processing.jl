using Pkg; Pkg.activate("..")

using OceananigansAnalysis, Glob, Printf

const noutput = 32

simname = "wind_stress"
rundir = "data"
curdir = pwd()

data_repo = "."

try
    cd(rundir)

    @info "Changed to $rundir"

    datafiles = glob("*.nc", ".")

    @assert length(datafiles) > 0 "No data files were found in $rundir"

    @info "The data files in here are \n $(("$file \n" for file in datafiles)...)"

    name = name_without_num(basename(datafiles[1]))
    timeseries_filepath = joinpath(name * "_timeseries.jld2")

    @assert !isfile(timeseries_filepath) "Found $timeseries_filepath in $rundir. Moving on..."

    @info "Creating a timeseries at $timeseries_filepath ..."

    try
        @time create_timeseries(timeseries_filepath, simname=simname, noutput=noutput, verbose=true)
    catch inner
        if :msg in propertynames(inner)
            @warn "Recieved $(typeof(inner)): $(inner.msg) \n Deleting $timeseries_filepath."
        else 
            @warn "Recieved $(typeof(inner))\n Deleting $timeseries_filepath."
        end

        rm(timeseries_filepath)
    end
catch outer
    if :msg in propertynames(inner)
        @warn "Recieved $(typeof(inner)): $(inner.msg) \n Deleting $timeseries_filepath."
    else 
        @warn "Recieved $(typeof(inner))\n Deleting $timeseries_filepath."
    end
finally
    cd(curdir)
end

