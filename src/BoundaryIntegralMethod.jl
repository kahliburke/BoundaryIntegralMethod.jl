__precompile__()
module BoundaryIntegralMethod
using Distributed

export Boundary, SimpleBoundary, Circle, DiscretizedBoundary, Point2D, nExt,
      nInt, isInside, nIdxRefr, boundaryPoint, boundaryPointDeriv, normVector, boundaryIdx,
      scaledParameter, absDiff, boundaryPointVector, A, A′, residualInverseIteration, NoeckelQuadrupoleBoundary,
      exportTotalField, ψField, G, ∂νG, arcLen, curvature, unitNormVector

include("boundaries.jl")
include("bimCalculations.jl")
include("modeSearch.jl")
include("totalField.jl")
include("DirectSweep.jl")
include("Follow.jl")
include("Husimi.jl")
include("Utilities.jl")

end # module
# @everywhere using BoundaryIntegralMethod
# @everywhere using DistributedArrays
