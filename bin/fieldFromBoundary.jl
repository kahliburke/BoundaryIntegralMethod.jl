#!/usr/bin/env julia
include("commandUtils.jl")
activatePkg()
using ArgParse
using JLD
using Distributed
s = ArgParseSettings(description="Calculates field based on boundary saved in file.")
@add_arg_table s begin
    "--simFile", "-s"
        help = "Simulation file to include, must set up boundary named simBoundary, can also edit or set k or other params to avoid the need to
                set in command arguments."
        required = true
    "--boundaryFile", "-b"
        help = "File with calculated kf and field xf"
        required = true
    "--dx"
        help = "Plot dx"
        arg_type = Float64
        default = 0.01
    "--dy"
        help = "Plot dy"
        arg_type = Float64
        default = 0.01
    "--xRange"
        help = "X range (tuple)"
        default = "(-1.5, 1.5)"
    "--yRange"
        help = "Y range (tuple)"
        default = "(-1.5, 1.5)"
    "--procs", "-p"
        help = "Number of processes to use in parallel"
        arg_type = Int
        default = convert(Int, Sys.CPU_THREADS/2)
end

arguments = parse_args(s)
@debug "Arguments" arguments

addprocs(arguments["procs"])

using BoundaryIntegralMethod

boundaryFile = arguments["boundaryFile"]
xRange = eval(Meta.parse(arguments["xRange"]))
yRange = eval(Meta.parse(arguments["yRange"]))
dx = arguments["dx"]
dy = arguments["dy"]

simFile = arguments["simFile"]
simName = split(basename(simFile), ".")[1]
if (!startswith(simFile, "/"))
  simFile = "$(pwd())/$simFile"
end
include(simFile)
b = simBoundary

(k, xf) = jldopen(boundaryFile, "r") do file
  (read(file, "k/full"), read(file, "boundary/full"))
end

n = convert(Int,size(xf,1)/4)

kForFile = string(real(k),"_",imag(k))
outFile = "output/field-$simName-$n-$kForFile.h5"
mkpath("output")
exportTotalField(b, xf, k; outFile=outFile, xRange=xRange, yRange=yRange, dx=dx, dy=dy)
