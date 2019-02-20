#!/usr/bin/env julia
include("commandUtils.jl")
activatePkg()
using ArgParse
using JLD
using Distributed
s = ArgParseSettings(description="Resonant Mode Direct sweep, implements algorithm 3 in Heider paper.")
@add_arg_table s begin
    "--simFile", "-s"
        help = "Simulation file to include, must set up boundary named simBoundary, can also edit or set k or other params to avoid the need to
                set in command arguments."
        required = true
    "-n"
        help = "Discretizaton N"
        arg_type = Int
        default = 50
    "--k0"
        help = "Search for resonances close to this k (can also set k0 in simFile)."
        arg_type = Complex{Float64}
    "--procs", "-p"
        help = "Number of processes to use in parallel"
        arg_type = Int
        default = convert(Int, Sys.CPU_THREADS/2)
end
arguments = parse_args(s)

@debug arguments

addprocs(arguments["procs"])
using BoundaryIntegralMethod
using BoundaryIntegralMethod.DirectSweep
n = arguments["n"]
k0 = arguments["k0"]

simFile = arguments["simFile"]
simName = split(basename(simFile), ".")[1]
if (!startswith(simFile, "/"))
  simFile = "$(pwd())/$simFile"
end
include(simFile)
b = simBoundary

(ks, xs) = directSweep(b, n, k0)
mkpath("output")
jldopen("output/sweep-$simName-$n-$(real(k0))_$(imag(k0))im.h5", "w") do file
  write(file, "ks/full", ks)
  write(file, "ks/Re", real(ks))
  write(file, "ks/Im", imag(ks))
  write(file, "xs/full", xs)
  write(file, "xs/Re", real(xs))
  write(file, "xs/Im", imag(xs))
end
exit(0)
