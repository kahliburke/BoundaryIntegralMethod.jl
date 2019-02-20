using JLD
using ProgressMeter
using Distributed

# See Wiersig paper for details/definitions of these functions. There is an overall sign change
# in the Green's function to match definitions in Heider's paper.
"""
    ∂νG(b::Boundary, r′::Point2D, t::Float64, n::Complex{Float64}, k::Complex{Float64})

Normal derivative of Green's function at point r′ (external or internal to boundary),
point parameterized by t on the boundary, and specified index of refraction and complex
wave vector k.
"""
function ∂νG(b::Boundary, r′::Point2D, t::Float64, n::Complex{Float64}, k::Complex{Float64})
  r = boundaryPoint(b, t)
  diff = r-r′
  absDiff = norm(diff)
  nv = normVector(b, t)
  # Normalize this normal vector because Wiersig paper defines it as unit vector.
  ν = nv/norm(nv)
  cosα = dot(Point2D(ν), diff)/absDiff
  -im*n*k/4*cosα*BoundaryIntegralMethod.h1(n*k*absDiff)
end

"""
    G(b::Boundary, r′::Point2D, t::Float64, n::Complex{Float64}, k::Complex{Float64})

Green's function at point r′ (external or internal to boundary),
point parameterized by t on the boundary, and specified index of refraction and complex
wave vector k.

"""
function G(b::Boundary, r′::Point2D, t::Float64, n::Complex{Float64}, k::Complex{Float64})
  r = boundaryPoint(b, t)
  absDiff = norm(r-r′)
  1/4*im*BoundaryIntegralMethod.h0(n*k*absDiff)
end

function ds(b::Boundary, j::Int, N::Int)
  arcLen(b, j*π/N, 2N)
end

"""
    ψField(b::Boundary, N::Int, r′::Point2D, n::Complex{Float64}, k::Complex{Float64}, ψ_Γ::Array{Complex{Float64},1}, ∂νψ_Γ::Array{Complex{Float64},1})

Calculate the field at a point given
the field and its normal derivative on the boundary.
"""
function ψField(b::Boundary, N::Int, ptSpacing::Float64, r′::Point2D, n::Complex{Float64}, k::Complex{Float64}, ψ_Γ::Array{Complex{Float64},1}, ∂νψ_Γ::Array{Complex{Float64},1}, dsList::Array{Float64})
  # Note that due to sign differences in the Green's function in Heider's paper, there is a plus in the discrete integral below
  # instead of the customary minus
  fieldVal = try
    sum([dsList[j+1]*(ψ_Γ[j+1]*∂νG(b, r′, j*π/N, n, k) + G(b, r′, j*π/N, n, k)*∂νψ_Γ[j+1]) for j in 0:2N-1])
  catch err
    @error exception=err
    foundPointIdx = -1
    for j in 0:2N-1
      pt = boundaryPoint(b, j*π/N)
      if isapprox(norm(pt-r′), 0.0, atol=ptSpacing/4)
        foundPointIdx = j
        break
      end
    end
    @info "Singular field at r′ = $r′, foundPointIdx = $foundPointIdx, t = $(foundPointIdx*π/N), boundary point: $(boundaryPoint(b, foundPointIdx*π/N)), using boundary value $(ψ_Γ[foundPointIdx+1])"
    ψ_Γ[foundPointIdx+1]
  end

  # fieldVal
  converted = Base.convert(Complex{Float32}, fieldVal)
end

function exportTotalField(b::Boundary, computedBoundary::Array{Complex{Float64},1}, k::Complex{Float64};
    outFile::String = nothing, xRange=(-1.5,1.5), yRange=(-1.5,1.5), dx=.01, dy=.01, boundaryPoints = 200)
  @info "Calculating field in range ($xRange, $yRange) dx=$dx, dy=$dy..."
  ys = yRange[1]:dy:yRange[2]
  xs = xRange[1]:dx:xRange[2]
  #field = Array(Complex{Float64},(length(ys),length(xs)))
  N = Base.convert(Int,length(computedBoundary)/4)
  ψ_Γ = computedBoundary[1:2N]
  ∂νψ_Γ = computedBoundary[2N+1:end]
  dr = sqrt(dx^2+dy^2)
  discreteB = DiscretizedBoundary(b, boundaryPoints)

  dsList = [ds(b, j, N) for j in 0:2N-1]
  fieldRows = []
  @showprogress 1 for i in 1:length(ys)
    row = @distributed (hcat) for j in 1:length(xs)
      r′ = Point2D(xs[j], ys[i])
      inside = isInside(discreteB, r′)
      n = inside ? nInt(b) : nExt(b)
      ψField(b, N, dr, r′, n, k, ψ_Γ, ∂νψ_Γ, dsList)
    end
    push!(fieldRows, row)
  end
  field = vcat(fieldRows...)

  if outFile != nothing
    JLD.jldopen(outFile, "w") do file
      write(file, "meta/k/full", k)
      write(file, "meta/k/Re", real(k))
      write(file, "meta/k/Im", imag(k))
      write(file, "meta/xRange", xRange)
      write(file, "meta/yRange", yRange)
      write(file, "meta/boundaryObject", b)
      write(file, "boundary/full", computedBoundary)
      write(file, "boundary/Re", real(computedBoundary))
      write(file, "boundary/Im", imag(computedBoundary))
      write(file, "field/full", field)
      write(file, "field/Re", real(field))
      write(file, "field/Im", imag(field))
    end
  end
  field
end
