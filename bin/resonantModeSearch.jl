#!/usr/bin/env julia
include("commandUtils.jl")
activatePkg()
using ArgParse
using JLD
using Distributed
s = ArgParseSettings(description="Resonant Mode Search, runs algorithms 1 & 2 from Heider paper.")
@add_arg_table s begin
    "--simFile", "-s"
        help = "Simulation file to include, must set up boundary named simBoundary, can also edit or set k or other params to avoid the need to
                set in command arguments."
        required = true
    "--fieldFile", "-f"
        help = "File with calculated k's and field approximation as initial guess for x0."
    "--fieldIdx", "-i"
        help = "Index in field file to use for k and x0."
        arg_type = Int
        default = 1
    "--k0"
        help = "Search for resonances close to this k (can also set k0 in simFile)."
        arg_type = Complex{Float64}
    "--eps"
        help = "Epsilon value, stop iteration when k changes less than this value "
        arg_type = Float64
        default = .000001
    "--dx"
        help = "Total field dx"
        arg_type = Float64
        default = 0.01
    "--dy"
        help = "Total field dy"
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

@debug arguments

addprocs(arguments["procs"])
@everywhere using BoundaryIntegralMethod

ϵ = arguments["eps"]
fieldFile = arguments["fieldFile"]
fieldIdx = arguments["fieldIdx"]

xRange = eval(Meta.parse(arguments["xRange"]))
yRange = eval(Meta.parse(arguments["yRange"]))
dx = arguments["dx"]
dy = arguments["dy"]

simFile = arguments["simFile"]
simName = importSim(simFile)
b = simBoundary

ks = nothing
xs = nothing
if (fieldFile != nothing)
  (ks, xs) = jldopen(fieldFile, "r") do file
    (read(file, "ks/full"), read(file, "xs/full"))
  end
end

if (ks != nothing)
  k0 = ks[fieldIdx]
  x0 = xs[fieldIdx,:]
end

n = convert(Int,size(x0,1)/4)

@info "Looking for resonance at k=$k0"

(kf, xf) = residualInverseIteration(10, ϵ, b, n, k0, x0)

@info "Final k=$kf"

kfForFile = string(real(kf),"_",imag(kf))
outFile = "output/field-$simName-$n-$kfForFile.h5"
mkpath("output")
boundaryOutFile = "output/boundaryField-$simName-$n-$kfForFile.h5"
jldopen(boundaryOutFile, "w") do file
  write(file, "k/full", kf)
  write(file, "k/Re", real(kf))
  write(file, "k/Im", imag(kf))
  write(file, "boundary/full", xf)
  write(file, "boundary/Re", real(xf))
  write(file, "boundary/Im", imag(xf))
end

exportTotalField(b, xf, kf; outFile=outFile, xRange=xRange, yRange=yRange, dx=dx, dy=dy)
