module Follow
using Plots
using BoundaryIntegralMethod
using BoundaryIntegralMethod.DirectSweep
using JLD

function followResonance(k0, b0, changeFn, steps; convergence = 1e-6, n=200, outputDir=nothing, fileBase="",
    maxDist=nothing, plotFn = nothing, startStep=1, namingFn = nothing, callbackFn = nothing)
  plots = []
  kfs = []
  bs = []
  b = b0
  k = k0
  if outputDir == nothing
    outputDir = "output/follow-$k0"
  end
  if fileBase != ""
    fileBase = fileBase+"-"
  end
  mkpath(outputDir)
  if plotFn == nothing
    plotFn = (b, xf, kf, baseName) -> exportTotalField(b, xf, kf; outFile="$baseName.h5", xRange=(-2,2), yRange=(-2,2), dx=.01, dy=.01);
  end

  for i in startStep:steps
    try
      (ks, xs) = directSweep(b, n, k)
      ks = ks[1:6]
      # Try closest resonance
      @debug ks
      p = sortperm(abs.(imag.(ks)-imag(k)))
      @info "Choosing from ks:"
      [@info "$i: $(ks[i])" for i in 1:length(ks)]
      # candidateK = ks[p[1]]
      # candidateX = xs[p[1],:]
      candidateK = ks[1]
      candidateX = xs[1,:]
      @info "Sortperm would be $(p[1]): $(ks[p[1]]), abs would be 1: $(ks[1])"
      @info "Choosing 1: $(ks[1])"

      (kf, xf) = residualInverseIteration(25, convergence, b, n, candidateK, candidateX)
      @info "At step $i, converged to $kf"
      paramName = if (namingFn != nothing)
        "-$(namingFn(b))"
      else
        ""
      end

      baseName= "$outputDir/$fileBase$i$paramName-$kf"
      field = plotFn(b, xf, kf, baseName)
      p = heatmap(clamp.(abs.(field), 0, .6), cbar=false, aspect_ratio=1, ticks=nothing, title="$i$paramName-k=$kf")
      savefig("$baseName.png")
      if maxDist != nothing && maxDist < abs(k-kf)
        @error "Stopping after $i steps, max distance between wavenumber exceeded in a step ($k -> $kf)"
        break
      end
      push!(kfs, kf)
      push!(plots, p)
      push!(bs, b)
      writeDatafile("$outputDir/$(fileBase)data.jld", kfs, bs)
      if callbackFn != nothing
        callbackFn(b, kf, xf, "$fileBase$i$paramName", baseName)
      end
      k = kf
      b = changeFn(b, i)
    catch e
      @error "Stopping after $i steps, caught error" exception=e
      break
    end
  end
  (kfs, plots)
end

function writeDatafile(name, kfs, bs)
  jldopen(name, "w") do file
    write(file, "ks/full", kfs)
    write(file, "ks/Re", real(kfs))
    write(file, "ks/Im", imag(kfs))
    write(file, "boundaries", bs)
  end
end

end
