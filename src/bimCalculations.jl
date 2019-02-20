using SpecialFunctions
using DistributedArrays
using LinearAlgebra

mutable struct BoundaryPrecomputedCalcs
  b::Boundary
  # Number of discretized steps (values go from 0:(2N-1))
  N::Int32
  k::Complex{Float64}
  j0NiMatrix::Array{Complex{Float64},2}
  j0NeMatrix::Array{Complex{Float64},2}
  j1NiMatrix::Array{Complex{Float64},2}
  j1NeMatrix::Array{Complex{Float64},2}
  h0NiMatrix::Array{Complex{Float64},2}
  h0NeMatrix::Array{Complex{Float64},2}
  h1NiMatrix::Array{Complex{Float64},2}
  h1NeMatrix::Array{Complex{Float64},2}
  absDiff::Array{Float64,2}
end

function absDiff(b::Boundary, t::Float64, τ::Float64)
  a = boundaryPoint(b, t)
  b = boundaryPoint(b, τ)
  norm(a-b)
end

function absDiff(b::Boundary, tτMatrix::Array{Array{Float64,1},2})
  map(pair -> norm(boundaryPoint(b, pair[1])-boundaryPoint(b, pair[2])), tτMatrix)
end

# @param l index (from 1:(2N-1) )
function γ(b::Boundary, N::Int64, l::Int64)
  norm(boundaryPointDeriv(b, l*π/N))+0im
end

# Quadrature weight factor
function R(j::Int64, N::Int64)
  cosAry = Array{Float64}(undef, N-1)
  for m in 1:N-1
    cosAry[m] = 1/m * cos(m*j*π/N)
  end
  -2π/N * sum(cosAry) - (-1)^j*π/N^2
end

################
j0(x::Complex{Float64}) = besselj(0, x)
j1(x::Complex{Float64}) = besselj(1, x)
h0(x::Complex{Float64}) = hankelh1(0, x)
# this is more accurate for some small values but it didn't seem to help me
# h0(x::Complex{Float64}) = j0(x)+bessely(0, x)*im
h1(x::Complex{Float64}) = hankelh1(1, x)
# h1(x::Complex{Float64}) = j1(x)+bessely(1, x)*im
nVals(b::Boundary, t::Float64) = (nInt(b), nExt(b))
xVals(b::Boundary, t::Float64, τ::Float64) = (convert(Array{Float64,1}, boundaryPoint(b, t)), convert(Array{Float64,1}, boundaryPoint(b, τ)))
log4Term(t::Float64, τ::Float64)::Float64 = log(4*sin((t-τ)/2)^2)

################
# H functions
function H1(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    # Note: Need to check that we can use the k values at t, not τ, this is a deviation from the Heider paper
    (ni, ne) = nVals(b, t)
    (xt, xτ) = xVals(b, t, τ)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j1ni = bpc.j1NiMatrix[tIdx, τIdx]
    j1ne = bpc.j1NeMatrix[tIdx, τIdx]
    -1/(2π)*dot(normVector(b, τ), (xt-xτ)/ad)*(k*ne*j1ne-k*ni*j1ni)
  end
end

function H2(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    H(b, t, τ, k, bpc)-H1(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function H(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (xt, xτ) = xVals(b, t, τ)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h1ni = bpc.h1NiMatrix[tIdx, τIdx]
  h1ne = bpc.h1NeMatrix[tIdx, τIdx]
  im/2*dot(normVector(b, τ), (xt-xτ)/ad)*(k*ne*h1ne-k*ni*h1ni)
end

H1star(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H1(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))
H2star(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H2(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))
Hstar(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))

###################
# L functions
function L1(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  if (t == τ)
    ((k*ni)^2-(k*ne)^2) / (2π)
  else
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    -1/(2π)*((k*ne)^2*j0ne-(k*ni)^2*j0ni)
  end
end

function L2(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    (ni, ne) = nVals(b, t)
    (im/2 - Base.MathConstants.eulergamma/π - 1/π*log(norm(boundaryPointDeriv(b, t)) / 2)) * ((k*ne)^2-(k*ni)^2) - 1/π*((k*ne)^2*log(k*ne) - (k*ni)^2*log(k*ni))
  else
    L(b, t, τ, k, bpc)-L1(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function L(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  im/2*((k*ne)^2*h0ne-(k*ni)^2*h0ni)
end

###################
# M functions
function M1(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    (ni, ne) = nVals(b, t)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    -1/(2π)*(j0ne-j0ni)
  end
end

function M2(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    (ni, ne) = nVals(b, t)
    1/π*log(ni/ne)
  else
    M(b, t, τ, k, bpc) - M1(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function M(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  im/2 * (h0ne-h0ni)
end

###################
# N functions
function N(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  xtP = boundaryPointDeriv(b, t)
  xτP = boundaryPointDeriv(b, τ)
  xpdot = dot(xtP, xτP)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  h1ni = bpc.h1NiMatrix[tIdx, τIdx]
  h1ne = bpc.h1NeMatrix[tIdx, τIdx]
  nbar = Nbar(b, t, τ, bpc)
  im/2*nbar*((k*ne)^2*h0ne-2*k*ne*h1ne/ad) +
    im/2*k*ne*xpdot/ad*h1ne - im/2*nbar*((k*ni)^2*h0ni-2*k*ni*h1ni/ad) -
    im/2*k*ni*xpdot/ad*h1ni
end

function Nbar(b::Boundary, t::Float64, τ::Float64, bpc::BoundaryPrecomputedCalcs)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  (xt, xτ) = xVals(b, t, τ)
  xtP = boundaryPointDeriv(b, t)
  xτP = boundaryPointDeriv(b, τ)

  xtm = (xt-xτ)
  dot(xtm, xtP)*dot(xtm, xτP)/ad^2
end

function N1(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  xtP = boundaryPointDeriv(b, t)

  if (t == τ)
    norm(xtP)^2*((k*ni)^2-(k*ne)^2)/(4π)
  else
    xτP = boundaryPointDeriv(b, τ)
    xpdot = dot(xtP, xτP)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    j1ni = bpc.j1NiMatrix[tIdx, τIdx]
    j1ne = bpc.j1NeMatrix[tIdx, τIdx]

    -1/(2π)*Nbar(b, t, τ, bpc)*((k*ne)^2*j0ne-2*k*ne*j1ne/ad) -
      k*ne*xpdot/(2π*ad)*j1ne +
      1/(2π)*Nbar(b, t, τ, bpc)*((k*ni)^2*j0ni-2*k*ni*j1ni/ad) +
      k*ni*xpdot/(2π*ad)*j1ni
  end
end

function N2(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    (xt, xτ) = xVals(b, t, τ)
    (ni, ne) = nVals(b, t)
    xtP = boundaryPointDeriv(b, t)
    normXtP = norm(xtP)
    normXtPSq = norm(xtP)^2

    normXtPSq/(4π) * (
      ((k*ne)^2-(k*ni)^2)*(π*im-1-2*Base.MathConstants.eulergamma) -
        2*(k*ne)^2*log(k*ne*normXtP/2)+2*(k*ni)^2*log(k*ni*normXtP/2)
    )
  else
    N(b, t, τ, k, bpc)-N1(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

###################
# A matrix blocks
#
A11(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ -(R(abs(j-l), N)*H1(b, j*π/N, l*π/N, k, bpc) + π/N*H2(b, j*π/N, l*π/N, k, bpc)) for j in 0:(2N-1), l in 0:(2N-1) ]

function A12(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  @DArray [ -(R(abs(j-l), N)*M1(b, j*π/N, l*π/N, k, bpc)*γ(b, N, l) +
    π/N*M2(b, j*π/N, l*π/N, k, bpc)*γ(b,N,l))
    for j in 0:(2N-1), l in 0:(2N-1) ]
end

function A21entry(b::Boundary, N::Int64, k::Complex{Float64}, j::Int64, l::Int64, bpc::BoundaryPrecomputedCalcs)
  xjP = boundaryPointDeriv(b, j*π/N)
  xlP = boundaryPointDeriv(b, l*π/N)
  R(abs(j-l), N) * (L1(b, j*π/N, l*π/N, k, bpc)*dot(xjP, xlP) - N1(b, j*π/N, l*π/N, k, bpc)) +
    π/N*(L2(b, j*π/N, l*π/N, k, bpc)*dot(xjP, xlP) - N2(b, j*π/N, l*π/N, k, bpc))
end
A21(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ A21entry(b, N, k, j, l, bpc) for j in 0:(2N-1), l in 0:(2N-1) ]

A22(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ R(abs(j-l), N)*H1star(b, j*π/N, l*π/N, k, bpc) + π/N*H2star(b, j*π/N, l*π/N, k, bpc) for j in 0:(2N-1), l in 0:(2N-1) ]

###################
# Other matrix blocks
γMatrix(b::Boundary, N::Int64) = Diagonal([γ(b, N, l) for l in 0:(2N-1)])

###################
# Full matrix block
function A(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  a11 = A11(b, N, k, bpc)
  a12 = A12(b, N, k, bpc)
  a21 = A21(b, N, k, bpc)
  a22 = A22(b, N, k, bpc)
  2*[ Matrix(1.0I, 2N, 2N) zeros(2N,2N); zeros(2N,2N) γMatrix(b, N)] - [convertAndClose(a11) convertAndClose(a12); convertAndClose(a21) convertAndClose(a22)]
end

function A(b::Boundary, N::Int64, k::Complex{Float64})
  bpc = doPrecalculation(b, N, k)
  A(b, N, k, bpc)
end

###################
# A'(k) derivative functions

################
# H functions
function H1′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    # Note: Need to check that we can use the k values at t, not τ, this is a deviation from the Heider paper
    (ni, ne) = nVals(b, t)
    (xt, xτ) = xVals(b, t, τ)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    k*ad/(2π)*dot(normVector(b, τ), (xt-xτ)/ad)*(ni^2*j0ni-ne^2*j0ne)
  end
end

function H2′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    H′(b, t, τ, k, bpc)-H1′(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function H′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (xt, xτ) = xVals(b, t, τ)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  im*k*ad/2*dot(normVector(b, τ), (xt-xτ)/ad)*(ne^2*h0ne-ni^2*h0ni)
end

H1star′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H1′(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))
H2star′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H2′(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))
Hstar′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = H′(b, τ, t, k, bpc)*norm(boundaryPointDeriv(b, τ))

###################
# L functions
function L1′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  if (t == τ)
    k*(ni^2-ne^2) / π
  else
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    j1ni = bpc.j1NiMatrix[tIdx, τIdx]
    j1ne = bpc.j1NeMatrix[tIdx, τIdx]
    1/(2π)*k*(k*ad*(ne^3*j1ne-ni^3*j1ni)-2*ne^2*j0ne+2*ni^2*j0ni)
  end
end

function L2′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    (ni, ne) = nVals(b, t)
    k/π*(-2*ne^2*log(k*ne)+2*ni^2*log(k*ni)-(ne-ni)*(ne+ni)*(2*log(norm(boundaryPointDeriv(b, t)) / 2)-π*im+2*Base.MathConstants.eulergamma+1))
  else
    L′(b, t, τ, k, bpc)-L1′(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function L′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  h1ni = bpc.h1NiMatrix[tIdx, τIdx]
  h1ne = bpc.h1NeMatrix[tIdx, τIdx]
  -k*im/2*(k*ad*(ne^3*h1ne-ni^3*h1ni)-2*ne^2*h0ne+2*ni^2*h0ni)
end

###################
# M functions
function M1′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    (ni, ne) = nVals(b, t)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j1ni = bpc.j1NiMatrix[tIdx, τIdx]
    j1ne = bpc.j1NeMatrix[tIdx, τIdx]
    ad/(2π)*(ne*j1ne-ni*j1ni)
  end
end

function M2′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    0.0+0im
  else
    M′(b, t, τ, k, bpc) - M1′(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

function M′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h1ni = bpc.h1NiMatrix[tIdx, τIdx]
  h1ne = bpc.h1NeMatrix[tIdx, τIdx]
  -im*ad/2 * (ne*h1ne-ni*h1ni)
end

###################
# N functions
function N′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  xtP = boundaryPointDeriv(b, t)
  xτP = boundaryPointDeriv(b, τ)
  xpdot = dot(xtP, xτP)
  (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
  ad = bpc.absDiff[tIdx, τIdx]
  h0ni = bpc.h0NiMatrix[tIdx, τIdx]
  h0ne = bpc.h0NeMatrix[tIdx, τIdx]
  h1ni = bpc.h1NiMatrix[tIdx, τIdx]
  h1ne = bpc.h1NeMatrix[tIdx, τIdx]
  k*im/2 * (
    k*ad*Nbar(b, t, τ, bpc) * (ni^3*h1ni-ne^3*h1ne) +
    xpdot*(ne^2*h0ne-ni^2*h0ni)
  )
end

function N1′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  (ni, ne) = nVals(b, t)
  xtP = boundaryPointDeriv(b, t)
  xτP = boundaryPointDeriv(b, τ)

  if (t == τ)
    norm(xtP)^2*k*(ni^2-ne^2)/(2π)
  else
    xτP = boundaryPointDeriv(b, τ)
    xpdot = dot(xtP, xτP)
    (tIdx, τIdx) = tAndτIndices(bpc, t, τ)
    ad = bpc.absDiff[tIdx, τIdx]
    j0ni = bpc.j0NiMatrix[tIdx, τIdx]
    j0ne = bpc.j0NeMatrix[tIdx, τIdx]
    j1ni = bpc.j1NiMatrix[tIdx, τIdx]
    j1ne = bpc.j1NeMatrix[tIdx, τIdx]
    k/(2π)*(
      k*ad*Nbar(b, t, τ, bpc)*(ne^3*j1ne-ni^3*j1ni) +
      xpdot*(ni^2*j0ni-ne^2*j0ne)
    )
  end
end

function N2′(b::Boundary, t::Float64, τ::Float64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  if (t == τ)
    (ni, ne) = nVals(b, t)
    xtP = boundaryPointDeriv(b, t)
    normXtP = norm(xtP)
    normXtPSq = norm(xtP)^2

    k*normXtPSq/(2π) * (
      -2*ne^2*log(k/2*ne*normXtP) +
      2*ni^2*log(k/2*ni*normXtP) -
      (2+2*Base.MathConstants.eulergamma-π*im)*(ne-ni)*(ne+ni)
    )
  else
    N′(b, t, τ, k, bpc)-N1′(b, t, τ, k, bpc)*log4Term(t, τ)
  end
end

###################
# A matrix blocks′
#
A11′(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ -(R(abs(j-l), N)*H1′(b, j*π/N, l*π/N, k, bpc) + π/N*H2′(b, j*π/N, l*π/N, k, bpc)) for j in 0:(2N-1), l in 0:(2N-1) ]

function A12′(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  @DArray [ -(R(abs(j-l), N)*M1′(b, j*π/N, l*π/N, k, bpc)*γ(b, N, l) +
    π/N*M2′(b, j*π/N, l*π/N, k, bpc)*γ(b,N,l))
    for j in 0:(2N-1), l in 0:(2N-1) ]
end

function A21entry′(b::Boundary, N::Int64, k::Complex{Float64}, j::Int64, l::Int64, bpc::BoundaryPrecomputedCalcs)
  xjP = boundaryPointDeriv(b, j*π/N)
  xlP = boundaryPointDeriv(b, l*π/N)
  R(abs(j-l), N) * (L1′(b, j*π/N, l*π/N, k, bpc)*dot(xjP, xlP) - N1′(b, j*π/N, l*π/N, k, bpc)) +
    π/N*(L2′(b, j*π/N, l*π/N, k, bpc)*dot(xjP, xlP) - N2′(b, j*π/N, l*π/N, k, bpc))
end
A21′(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ A21entry′(b, N, k, j, l, bpc) for j in 0:(2N-1), l in 0:(2N-1) ]

A22′(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs) = @DArray [ R(abs(j-l), N)*H1star′(b, j*π/N, l*π/N, k, bpc) + π/N*H2star′(b, j*π/N, l*π/N, k, bpc) for j in 0:(2N-1), l in 0:(2N-1) ]

###################
# Full matrix block
function A′(b::Boundary, N::Int64, k::Complex{Float64}, bpc::BoundaryPrecomputedCalcs)
  a11′ = A11′(b, N, k, bpc)
  a12′ = A12′(b, N, k, bpc)
  a21′ = A21′(b, N, k, bpc)
  a22′ = A22′(b, N, k, bpc)
  - [convertAndClose(a11′) convertAndClose(a12′); convertAndClose(a21′) convertAndClose(a22′)]
end

function A′(b::Boundary, N::Int64, k::Complex{Float64})
  bpc = doPrecalculation(b, N, k)
  A′(b, N, k, bpc)
end

###################
# Utility functions
function tAndτIndices(bpc::BoundaryPrecomputedCalcs, t::Float64, τ::Float64)
  (Base.round(Int, t*bpc.N/π)+1, Base.round(Int, τ*bpc.N/π)+1)
end

function convertAndClose(a::DistributedArrays.DArray)
  tmp = Base.convert(Array, a)
  close(a)
  tmp
end

function doPrecalculation(b::Boundary, N::Int64, k::Complex{Float64})
  range = 0:(2N-1)
  nis = [nVals(b, l*π/N)[1] for l in range]
  nes = [nVals(b, l*π/N)[2] for l in range]

  tτMatrix = [[l*π/N,j*π/N] for l in range, j in range]
  absDiffM = absDiff(b, tτMatrix)

  localNiArgs = k*nis.*absDiffM
  localNeArgs = k*nes.*absDiffM
  niArgs = distribute(localNiArgs)
  neArgs = distribute(localNeArgs)
  j0NiM = besselj.(0, niArgs)
  j0NeM = besselj.(0, neArgs)
  j1NiM = besselj.(1, niArgs)
  j1NeM = besselj.(1, neArgs)

  # Set diagonal to NaN to prevent exception for hankel functions
  [localNiArgs[i, i] = NaN for i in 1:2N]
  [localNeArgs[i, i] = NaN for i in 1:2N]
  close(niArgs)
  close(neArgs)
  niArgs = distribute(localNiArgs)
  neArgs = distribute(localNeArgs)
  h0NiM = h0.(niArgs)
  h0NeM = h0.(neArgs)
  h1NiM = h1.(niArgs)
  h1NeM = h1.(neArgs)
  close(niArgs)
  close(neArgs)
  bpc = BoundaryPrecomputedCalcs(b, N, k, convertAndClose(j0NiM), convertAndClose(j0NeM), convertAndClose(j1NiM), convertAndClose(j1NeM),
      convertAndClose(h0NiM), convertAndClose(h0NeM), convertAndClose(h1NiM), convertAndClose(h1NeM), absDiffM)
  # @everywhere gc()
  bpc
end
