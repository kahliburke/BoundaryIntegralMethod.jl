#!/usr/bin/env julia
include("commandUtils.jl")
activatePkg()
using ArgParse
using JLD
using Distributed
using Plots
s = ArgParseSettings(description="Plot Field From H5 Ouput, plots field from output of BIM field calculation.")
@add_arg_table s begin
    "--fieldFile", "-f"
        help = "File with calculated k's and field approximation as initial guess for x0."
    "--clamp", "-c"
        help = "Clamp max field value below this number (to reduce dynamic range)."
        arg_type = Float64
        default = 1.0
end
arguments = parse_args(s)
@debug arguments

fieldFile = arguments["fieldFile"]
clampMax = arguments["clamp"]
field = jldopen(fieldFile) do file
   read(file, "field/full")
end
gr()
@info "Max field value (abs) = $(maximum(abs.(field)))"
Plots.heatmap(clamp.(abs.(field), 0, clampMax),size=(1000,1000), cbar=false, aspect_ratio=1)
Plots.savefig("$fieldFile.png")
