abstract type Boundary end

abstract type SimpleBoundary <: Boundary end
# Assumption of SimpleBoundary is that it has fields
# ni::Complex{Float64} - Internal index of refraction
# ne::Complex{Float64} - External index of refraction

nInt(s::SimpleBoundary) = s.ni
nExt(s::SimpleBoundary) = s.ne

struct Point2D
  x::Float64
  y::Float64
end

import Base:-
import Base:+
import Base.*
import LinearAlgebra.dot
import LinearAlgebra.norm
convert(::Type{Array{Float64,1}}, p::Point2D) = [p.x, p.y]
function convert(::Type{Point2D}, a::Array{Float64,1})
  if length(a) != 2
    error("Only 2 element arrays can be converted to points")
  end
  Point2D(a[1],a[2])
end
function -(a::Point2D, b::Point2D)
  Point2D(convert(Array{Float64,1}, a)-convert(Array{Float64,1}, b))
end
function +(a::Point2D, b::Point2D)
  Point2D(convert(Array{Float64,1}, a)+convert(Array{Float64,1}, b))
end
function +(a::Point2D, b::Array{Float64,1})
  Point2D(convert(Array{Float64,1}, a)+b)
end
function *(a::Number, p::BoundaryIntegralMethod.Point2D)
  Point2D(a*p.x, a*p.y)
end
Point2D(a::Array{Float64}) = convert(Point2D, a)
function dot(a::Point2D, b::Point2D)
  dot(convert(Array{Float64,1}, a), convert(Array{Float64,1}, b))
end
function dot(a::Array{Float64,1}, b::Point2D)
  dot(a, convert(Array{Float64,1}, b))
end
function norm(a::Point2D)
  sqrt(dot(a, a))
end

"""
Return boundary point as array.
"""
function boundaryPointVector(b::Boundary, t::Float64)
  convert(Array{Float64, 1}, boundaryPoint(b, t))
end

"""
Return boundary normal vector as array
"""
function normVector(b::Boundary, t::Float64)
  xp = boundaryPointDeriv(b, t)
  [xp[2], -xp[1]]
end

"""
Return unit scaled boundary normal vector as array
"""
function unitNormVector(b::Boundary, t::Float64)
  n = normVector(b, t)
  n/norm(n)
end

# This doesn't work with some boundaries like criss cross, probably due to discontinuities
# in boundary derivative
# function arcLen(b::Boundary, t::Float64, steps::Int)
#   function ds(t′::Float64)
#     (dxdt, dydt) = boundaryPointDeriv(b, t′)
#     sqrt(dxdt^2 + dydt^2)
#   end
#   quadgk(ds, t-pi/steps, t+pi/steps)[1]
# end

function arcLen(b::Boundary, t::Float64, steps::Int)
  p1 = boundaryPoint(b, t)
  nextT = t+2π/steps
  if nextT > 2π
    nextT = nextT - 2π
  end
  p2 = boundaryPoint(b, nextT)
  norm(p2-p1)
end

function curvature(b::Boundary, t::Float64)
  xp = boundaryPointDeriv(b, t)
  xpp = boundaryPointSecondDeriv(b, t)
  (xp[1]xpp[2]-xp[2]xpp[1])/(xp[1]^2+xp[2]^2)^(3/2)
end

function nIdxRefr(b::Boundary, p::Point2D)
  isInside(b, p) ? nInt(b) : nExt(b)
end

struct Circle <: SimpleBoundary
  # Center of circle
  center::Point2D
  # Radius
  r::Float64
  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}
end

function Circle(r, ni, ne)
  Circle(Point2D(0,0), r, ni, ne)
end

function boundaryPoint(c::Circle, t::Float64)
  Point2D(c.center.x+cos(t)*c.r, c.center.y+sin(t)*c.r)
end

function boundaryPointDeriv(c::Circle, t::Float64)
  c.r*[-sin(t), cos(t)]
end

function boundaryPointSecondDeriv(c::Circle, t::Float64)
  -c.r*[cos(t), sin(t)]
end

# Inverse of radius of curvature: ρ = 1/κ
function curvature(c::Circle, t::Float64)
  1/c.r
end

function arcLen(c::Circle, t::Float64, steps::Int)
  c.r*2π/steps
end

function isInside(c::Circle, p::Point2D)
  p = p-c.center
  norm(p) < c.r
end

function boundaryField(c::Circle, N::Int, m::Int, k::Complex{Float64})
  xbar = [exp(im*m*j*π/N)*besselj(m, c.ni*k*c.r) for j in 0:2N-1]
  append!(xbar, [exp(im*m*j*π/N)*1/2*c.ni*k*(besselj(m-1,c.ni*k*c.r)-besselj(m+1,c.ni*k*c.r)) for j in 0:2N-1])
  xbar
end

struct DiscretizedBoundary <: Boundary
  underlying::Boundary
  N::Int64
  boundaryPointVectors::Array{Array{Float64}}
  boundaryPointDerivs::Array{Array{Float64}}
  normVectors::Array{Array{Float64}}
  absDiffMatrix::Array{Float64,2}
end

function DiscretizedBoundary(b::Boundary, N::Int64)
  boundaryPointVectors = [boundaryPointVector(b, j*π/N) for j in 0:(2N-1)]
  boundaryPointDerivs = [boundaryPointDeriv(b, j*π/N) for j in 0:(2N-1)]
  normVectors = [normVector(b, j*π/N) for j in 0:(2N-1)]
  absDiffMatrix = [norm(boundaryPoint(b, j*π/N)-boundaryPoint(b, l*π/N)) for j in 0:(2N-1), l in 0:(2N-1)]
  DiscretizedBoundary(b, N, boundaryPointVectors, boundaryPointDerivs, normVectors, absDiffMatrix)
end

function tToIdx(d::DiscretizedBoundary, t::Float64)
  Int64(round(d.N*t/π))+1
end

function boundaryPoint(d::DiscretizedBoundary, t::Float64)
  i = tToIdx(d, t)
  x = d.boundaryPointVectors[i]
  Point2D(x[1], x[2])
end

function boundaryPointVector(d::DiscretizedBoundary, t::Float64)
  i = tToIdx(d, t)
  d.boundaryPointVectors[i]
end

function boundaryPointDeriv(d::DiscretizedBoundary, t::Float64)
  i = tToIdx(d, t)
  d.boundaryPointDerivs[i]
end

function normVector(d::DiscretizedBoundary, t::Float64)
  i = tToIdx(d, t)
  d.normVectors[i]
end

function nInt(d::DiscretizedBoundary)
  nInt(d.underlying)
end

function nExt(d::DiscretizedBoundary)
  nExt(d.underlying)
end

function absDiff(d::DiscretizedBoundary, t::Float64, τ::Float64)
  i = tToIdx(d, t)
  j = tToIdx(d, τ)
  d.absDiffMatrix[i,j]
end

function isInside(b::DiscretizedBoundary, p::Point2D)
  try
    polygonWindingNumber(b.boundaryPointVectors, p) != 0
  catch e
    # Consider points on edge or vertex to be outside
    false
  end
end

function polygonWindingNumber(boundaryVecs::Array{Array{Float64}}, p::Point2D)
  # See: http://www.sciencedirect.com/science/article/pii/S0925772101000128
  # "The point in polygon problem for arbitrary polygons"
  # An implementation of Hormann-Agathos (2001) Point in Polygon algorithm
  r = BoundaryIntegralMethod.convert(Array{Float64, 1}, p)
  ω = 0
  @inline detq(q1,q2) = (q1[1]-r[1])*(q2[2]-r[2])-(q2[1]-r[1])*(q1[2]-r[2])
  @inline crossing(q1, q2) = (q1[2] < r[2]) != (q2[2] < r[2])
  @inline right_crossing(q1, q2) = (detq(q1,q2) > 0) == (q2[2] > q1[2])
  @inline modify_ω(q1, q2) = ω = ω + 2*(q2[2] > q1[2]) - 1

  if r == boundaryVecs[1]
    throw(VertexException())
  end
  for idx in 1:length(boundaryVecs)
    q1 = boundaryVecs[idx]
    q2 = idx == length(boundaryVecs) ? boundaryVecs[1] : boundaryVecs[idx+1]
    if q2[2] == r[2]
      if q2[1] == r[1]
          throw(VertexException())
      elseif (q1[2] == r[2]) && ((q2[1] > r[1]) == (q1[1] < r[1]))
          throw(EdgeException())
      end
    end
    if crossing(q1, q2)
      if q1[1] >= r[1]
        if q2[1] > r[1]
          modify_ω(q1, q2)
        elseif right_crossing(q1, q2)
          modify_ω(q1, q2)
        end
      elseif q2[1] > r[1] && right_crossing(q1, q2)
        modify_ω(q1, q2)
      end
    end
  end
  return ω
end

struct HeiderConcaveBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D
  # Parameters for (x(t), y(t)) = (a cos t + b cos 2t - c, d sin t)
  # Heider paper uses a=1, b=.65, c=-.65, d=1.5
  a::Float64
  b::Float64
  c::Float64
  d::Float64
  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}
end

function HeiderConcaveBoundary(ni::Complex{Float64}, ne::Complex{Float64})
  HeiderConcaveBoundary(Point2D(0,0), 1, .65, -.65, 1.5, ni, ne)
end

function shapeFunc(b::HeiderConcaveBoundary, t::Float64)
  Point2D(b.a*cos(t) + b.b*cos(2t) - b.c, b.d*sin(t))
end

function boundaryPoint(b::HeiderConcaveBoundary, t::Float64)
  b.center+shapeFunc(b, t)
end

function boundaryPointDeriv(b::HeiderConcaveBoundary, t::Float64)
  [-b.a*sin(t)-2b.b*sin(2t), b.d*cos(t)]
end

function isInside(b::HeiderConcaveBoundary, p::Point2D)
  # figure out what t the angle of the poing would be at
  p = p - b.center
  ϕ = atan(p.y, p.x)
  # Range is -π < ϕ < π, place in range of t
  t = ϕ > 0 ? ϕ : 2π+ϕ
  p = p-c.center
  norm(p) < norm(boundaryPoint(b, t))
end

struct BowtieBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D

  # Parameters for shape in polar coords r(ϕ) = r0(1+ϵ cos(2ϕ))
  r0::Float64
  ϵ::Float64

  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}
end

function BowtieBoundary(r0::Float64, ϵ::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  BowtieBoundary(Point2D(0,0), r0, ϵ, ni, ne)
end

function radius(b::BowtieBoundary, t::Float64)
  b.r0*(1+b.ϵ*cos(2t))
end

function boundaryPoint(b::BowtieBoundary, t::Float64)
  r = radius(b, t)
  b.center+Point2D([r*cos(t), r*sin(t)])
end

function boundaryPointDeriv(b::BowtieBoundary, t::Float64)
  [-(1+2*b.ϵ+3*b.ϵ*cos(2t))*sin(t),
  (1-2*b.ϵ+3*b.ϵ*cos(2t))*cos(t)]
end

function isInside(b::BowtieBoundary, p::Point2D)
  # figure out what t the angle of the Point2D would be at
  p = p - b.center
  ϕ = atan(p.y,p.x)
  # Range is -π < ϕ < π, place in range of t
  t = ϕ > 0 ? ϕ : 2π+ϕ
  norm(p) < norm(boundaryPoint(b, t)-b.center)
end

# Another quadrupole boundary with normalized area, used in publications by J. Noeckel
struct NoeckelQuadrupoleBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D

  # Parameters for shape in polar coords r(ϕ) = (1 + ϵ Cos[ϕ]) (1+ϵ^2 /2 )^(-1/2)
  ϵ::Float64

  # Normally this is 2 for the standard quadrupole boundary but it can be increased to
  # create more oscillations in the boundary
  cosFactor::Int

  # Scaling applied to x-axis lobes only
  xScale::Float64

  #Phase for angle argument, allows rotation of the shape
  phase::Float64

  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}

  # Scale up shape by this factor
  scale::Float64
end

function NoeckelQuadrupoleBoundary(ϵ::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  NoeckelQuadrupoleBoundary(Point2D(0,0), ϵ, 2, 1.0, 1.0, 0., ni, ne, 1.0)
end

function NoeckelQuadrupoleBoundary(ϵ::Float64, ni::Complex{Float64}, ne::Complex{Float64}, scale::Float64)
  NoeckelQuadrupoleBoundary(Point2D(0,0), ϵ, 2, 1.0, 0., ni, ne, scale)
end

function NoeckelQuadrupoleBoundary(ϵ::Float64, cosFactor::Int, phase::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  NoeckelQuadrupoleBoundary(Point2D(0,0), ϵ, cosFactor, 1.0, phase, ni, ne, 1.0)
end

function NoeckelQuadrupoleBoundary(center::Point2D, ϵ::Float64, cosFactor::Int, phase::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  NoeckelQuadrupoleBoundary(center, ϵ, cosFactor, 1.0, phase, ni, ne, 1.0)
end


function radius(b::NoeckelQuadrupoleBoundary, t::Float64)
  (1 + b.ϵ*cos(b.cosFactor*t+b.phase))*(1 + b.ϵ^2/2)^(-1/2)
end

function boundaryPoint(b::NoeckelQuadrupoleBoundary, t::Float64)
  r = radius(b, t)*b.scale
  b.center+Point2D([b.xScale*r*cos(t), r*sin(t)])
end

function boundaryPointDeriv(b::NoeckelQuadrupoleBoundary, t::Float64)
  sq2 = sqrt(2)

  sint = sin(t)
  cost = cos(t)

  b.scale*[
    -b.xScale*sq2*(sint + b.ϵ*cos(b.cosFactor*t + b.phase)*sint + b.cosFactor*b.ϵ*cost*sin(b.cosFactor*t + b.phase)),
    sq2*(cost + b.ϵ*cos(b.cosFactor*t + b.phase)*cost - b.cosFactor*b.ϵ*sint*sin(b.cosFactor*t + b.phase))
  ]/sqrt(2+b.ϵ^2)
end

# Another shape based on octopole, like a 4-leaf clover
struct CloverBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D

  # Deformation parameter
  ϵ::Float64

  # scaling applied to x lobes only
  xScale::Float64

  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}

  # Scale up shape by this factor
  scale::Float64
end

function CloverBoundary(ϵ::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  CloverBoundary(Point2D(0,0), ϵ, 9.0/8.0, ni, ne, 1.0)
end

function radius(b::CloverBoundary, t::Float64)
  1 + 2*b.ϵ*atan(3*(1+cos(4t)))
end

function boundaryPoint(b::CloverBoundary, t::Float64)
  r = radius(b, t)*b.scale
  b.center+Point2D([b.xScale*r*cos(t), r*sin(t)])
end

function boundaryPointDeriv(b::CloverBoundary, t::Float64)

  sint = sin(t)
  cost = cos(t)
  term1 = (1+2*b.ϵ*atan(6*cos(2t)^2))
  term2 = -24*b.ϵ*sin(4t)/(1+36*cos(2t)^4)

  b.scale*[
    -b.xScale*term1*sint-b.xScale*term2*cost,
    term1*cost-term2*sint
  ]
end

# Another shape based on octopole, like a 4-leaf clover
struct SkewedOctupole <: SimpleBoundary
  # Center of structure
  center::Point2D

  # Deformation amount
  ϵ::Float64

  # Deformation applied to one corner
  η::Float64

  # t0 sets placement of deformation
  t0::Float64

  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}

  # Scale up shape by this factor
  scale::Float64
end

function SkewedOctupole(ϵ::Float64, η::Float64, t0::Float64, ni::Complex{Float64}, ne::Complex{Float64})
  SkewedOctupole(Point2D(0,0), ϵ, η, t0, ni, ne, 1.0)
end

function radius(b::SkewedOctupole, t::Float64)
  1 + b.ϵ*(1-b.η*((1+cos(t-b.t0))/2)^4)*cos(4t)
end

function boundaryPoint(b::SkewedOctupole, t::Float64)
  r = radius(b, t)*b.scale
  b.center+Point2D([r*cos(t), r*sin(t)])
end

function boundaryPointDeriv(b::SkewedOctupole, t::Float64)

  sintt0 = sin(t-b.t0)
  costt0 = cos(t-b.t0)
  cos4t = cos(4t)
  sin4t = sin(4t)
  term1 = 1/4*b.ϵ*b.η*cos4t*sintt0*(costt0+1)^3 -
    4b.ϵ*sin4t*(1-1/16*b.η*(costt0+1)^4)
  term2 = b.ϵ*cos4t*(1-1/16*b.η*(costt0+1)^4)+1

  b.scale*[
    cos(t)*term1-sin(t)*term2
    cos(t)*term2+sin(t)*term1
  ]
end

# Needed for Wiersig
# function boundaryPointSecondDeriv(b::NoeckelQuadrupoleBoundary, t::Float64)
#   sq2 = sqrt(2)
#
#   sint = sin(t)
#   cost = cos(t)
#
#   -b.scale*sq2/sqrt(2+b.ϵ^2)*[
#     cost + (1+b.cosFactor^2)*b.ϵ*cos(b.cosFactor*t + b.phase)*cost - 2*b.cosFactor*b.ϵ*sint*sin(b.cosFactor*t + b.phase),
#     sint + (1+b.cosFactor^2)*b.ϵ*cos(b.cosFactor*t + b.phase)*sint + 2*b.cosFactor*b.ϵ*cost*sin(b.cosFactor*t + b.phase)
#   ]
#
# end

# Another quadrupole boundary with normalized area, used in publications by J. Noeckel
struct ReuleauxBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D
  # Scale radius
  scale::Float64
  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}
  # (n, a_n) tuple coefficients defining the parameterized boundary shape
  coefficientTuples::Array{Tuple{Int, Complex},1}
end

function B1ReuleauxBoundary(ni::Complex{Float64}, ne::Complex{Float64}; scale::Float64=1.0/12)
  # B1 coefficients from 1539-3755/2014/90(2)/022903(15)
  b1Coeffs = [12., 0., 0., 3.0im/2, 0., 3.]
  rest = b1Coeffs[2:end]
  coeffTuples = vcat((0, b1Coeffs[1]), [(-n, conj(an)) for (n, an) in enumerate(rest)], [(n,an) for (n, an) in enumerate(rest)])
  ReuleauxBoundary(Point2D(0.,0.), scale, ni, ne, coeffTuples)
end

# The complex parameter describing the shape, where (x, y) = (re(z), im(z))
function z(r::ReuleauxBoundary, t::Float64)
  # Check for an != 0 is so we don't divide by 0 for a_1
  -im*r.coefficientTuples[1][2]-im*sum([an/(n+1)*(exp(im*(n+1)*t)-1) for (n, an) in r.coefficientTuples if an != 0])
end

# dz/dt
function z′(r::ReuleauxBoundary, t::Float64)
  # Check for an != 0 is so we don't divide by 0 for a_1
  sum([an*exp(im*(n+1)*t) for (n, an) in r.coefficientTuples if an != 0])
end

function boundaryPoint(r::ReuleauxBoundary, t::Float64)
  zPt = z(r, t)
  r.center+r.scale*Point2D([real(zPt), imag(zPt)])
end

function boundaryPointDeriv(r::ReuleauxBoundary, t::Float64)
  z′Pt = z′(r,t)
  r.scale*[real(z′Pt), imag(z′Pt)]
end

struct CrissCrossBoundary <: SimpleBoundary
  # Center of structure
  center::Point2D

  # Deformation parameter
  ϵ::Float64

  # Internal index of refraction
  ni::Complex{Float64}
  # External index of refraction
  ne::Complex{Float64}

  # Scale up shape by this factor
  scale::Float64

  # Length of right arm
  rLen::Float64

  # Length of left arm
  lLen::Float64

  tFlat::Float64

  # Shift parameter makes arm thinner but leaves the lobeSpeed
  shift::Float64
  CrissCrossBoundary(center, ϵ, ni, ne, scale, rLen, lLen) = new(center, ϵ, ni, ne, scale, rLen, lLen, tFlat(ϵ), 0.0)
  CrissCrossBoundary(center, ϵ, ni, ne, scale, rLen, lLen, shift) = new(center, ϵ, ni, ne, scale, rLen, lLen, tFlat(ϵ), shift)
  CrissCrossBoundary(ϵ, ni, ne) = new(Point2D(0,0), ϵ, ni, ne, 1.0, 2.0, 1.7, tFlat(ϵ), 0.0)
end

function rVec(b::CrissCrossBoundary, t::Float64)
  [cos(t), sin(t)]*(1+b.ϵ*cos(4t))
end

# Parameter at which shape becomes horizontal
function tFlat(ϵ::Float64)
  numeratorFactor = (sqrt(2)*sqrt(ϵ*(-5+13*ϵ)))/ϵ
  atan(sqrt(6+numeratorFactor)/(2*sqrt(5)),
    sqrt(14-numeratorFactor)/(2*sqrt(5)))
end

# Define squircle shape
function squircle(a::Float64, b::Float64, δ::Float64, sign1::Int,
    sign2::Int, t::Float64)
  [sign1*a*abs((sign1*cos(t)))^(1/2)+δ, sign2*b*abs(sin(sign2*t))^(1/2)]
end

function rEllipse(a::Float64, b::Float64, t::Float64)
  a*b/sqrt(b^2*cos(t)^2+a^2*sin(t)^2)
end

function circleMap(a::Float64, b::Float64, t::Float64)
  π/2*(a-rEllipse(a, b, t))/(a-b)
end

function radius(b::CrissCrossBoundary, t::Float64)
  (x0, y0) = rVec(b, b.tFlat)
  if b.tFlat <= t <= (π-b.tFlat)
    [0.0, -b.shift/2]+rVec(b, t)
  elseif (π+b.tFlat) <= t <= (2π-b.tFlat)
    [0.0, b.shift/2]+rVec(b, t)
  elseif 0 <= t <= b.tFlat
    squircle(b.rLen, y0-b.shift/2, x0, 1, 1, circleMap(b.rLen, y0, t/b.tFlat*π/2))
  elseif (π-b.tFlat) <= t <= π
    squircle(b.lLen, y0-b.shift/2, -x0, -1, 1, π-circleMap(b.lLen, y0, (π-t)/b.tFlat*π/2))
  elseif π <= t <= π+b.tFlat
    squircle(b.lLen, y0-b.shift/2, -x0, -1, -1, π+circleMap(b.lLen, y0, (t-π)/b.tFlat*π/2))
  else
    squircle(b.rLen, y0-b.shift/2, x0, 1, -1, 2π-circleMap(b.rLen, y0, (2π-t)/b.tFlat*π/2))
  end
end

function boundaryPoint(b::CrissCrossBoundary, t::Float64)
  r = radius(b, t)*b.scale
  b.center+Point2D([r[1], r[2]])
end

function lobeTangent(b::CrissCrossBoundary, t::Float64)
  denom = sqrt(4+34b.ϵ^2+8b.ϵ*cos(4t)-30b.ϵ^2*cos(8t))
  [(-2*sin(t)-3*b.ϵ*sin(3t)-5*b.ϵ*sin(5t))/denom,
    (2*cos(t)-3*b.ϵ*cos(3t)+5*b.ϵ*cos(5t))/denom]
end

function squircleTangent(a::Float64, b::Float64, δ::Float64, sign1::Int,
    sign2::Int, u::Float64)
  # Avoid issue with negative 0 sometimes being passed in which is outside domain
  u = abs(u)
  try
    if sign1 == 1 && sign2 == 1
      [-a/sqrt(abs(a^2+b^2*cot(u)^3)), b/sqrt(abs(b^2+a^2*tan(u)^3))]
    elseif sign1 == -1 && sign2 == 1
      [-a/sqrt(abs(a^2-b^2*cot(u)^3)), -b/sqrt(abs(b^2-a^2*tan(u)^3))]
    elseif sign1 == -1 && sign2 == -1
      [a/sqrt(abs(a^2+b^2*cot(u)^3)), -b/sqrt(abs(b^2+a^2*tan(u)^3))]
    else
      [a/sqrt(abs(a^2-b^2*cot(u)^3)), b/sqrt(abs(b^2-a^2*tan(u)^3))]
    end
  catch e
    @error "Domain error? $a, $b, $δ, $sign1, $sign2, $u\n$e" exception=e
  end
end

# Normalized tangent, needs speed factor to work for Heider
function normTangent(b::CrissCrossBoundary, t::Float64)
  (x0, y0) = rVec(b, b.tFlat)
  if b.tFlat <= t < (π-b.tFlat) || (π+b.tFlat) <= t < (2π-b.tFlat)
    lobeTangent(b, t)
  elseif 0 <= t < b.tFlat
    squircleTangent(b.rLen, y0-b.shift/2, x0, 1, 1, circleMap(b.rLen, y0, t/b.tFlat*π/2))
  elseif (π-b.tFlat) <= t < π
    squircleTangent(b.lLen, y0-b.shift/2, -x0, -1, 1, π-circleMap(b.lLen, y0, (π-t)/b.tFlat*π/2))
  elseif π <= t < π+b.tFlat
    squircleTangent(b.lLen, y0-b.shift/2, -x0, -1, -1, π+circleMap(b.lLen, y0, (t-π)/b.tFlat*π/2))
  else
    squircleTangent(b.rLen, y0-b.shift/2, x0, 1, -1, 2π-circleMap(b.rLen, y0, (2π-t)/b.tFlat*π/2))
  end
end

function squircleSpeed(a::Float64, b::Float64, t::Float64, tFlat::Float64)
  prefactor = 1/tFlat*π/2
  trigFactor = π*(a - a*b/sqrt(b^2*cos(t)^2+a^2*sin(t)^2))/(2*(a-b))
  spd = if 0 < t < π/2 && !isapprox(0, trigFactor) && !isapprox(π/2, trigFactor)
    spdSquared = (a+b)^2*π^2*cos(t)^2*sin(t)^2*
      (a^2*b^4*cos(trigFactor)*cot(trigFactor)+a^4*b^2*sin(trigFactor))*tan(trigFactor) /
      (16*(b^2*cos(t)^2 + a^2*sin(t)^2)^3)
      sqrt(spdSquared)
  elseif t <= 0 || !isapprox(0, trigFactor)
    sqrt(1/4*a*(a+b)*π)
  else
    sqrt(1/4*b*(a+b)*π)
  end
  prefactor*spd
end

function lobeSpeed(b::CrissCrossBoundary, t::Float64)
  sqrt((2+17*b.ϵ^2+4*b.ϵ*cos(4t)-15*b.ϵ^2*cos(8t))/2)
end

function boundaryPointDeriv(b::CrissCrossBoundary, t::Float64)
  (x0, y0) = rVec(b, b.tFlat)
  speed = if b.tFlat <= t <= (π-b.tFlat) || (π+b.tFlat) <= t <= (2π-b.tFlat)
    lobeSpeed(b, t)
  elseif 0 <= t <= b.tFlat
    squircleSpeed(b.rLen, y0-b.shift/2, t/b.tFlat*π/2, b.tFlat)
  elseif (π-b.tFlat) <= t <= π
    squircleSpeed(b.lLen, y0-b.shift/2, (π-t)/b.tFlat*π/2, b.tFlat)
  elseif π <= t <= π+b.tFlat
    squircleSpeed(b.lLen, y0-b.shift/2, (t-π)/b.tFlat*π/2, b.tFlat)
  else
    squircleSpeed(b.rLen, y0-b.shift/2, (2π-t)/b.tFlat*π/2, b.tFlat)
  end
  speed*normTangent(b, t)
end
