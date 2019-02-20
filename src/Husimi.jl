module Husimi
export husimiDistribution, husimi, coherentState
using BoundaryIntegralMethod
using DistributedArrays

"""
Represents the coherent state along the boundary parameterized by

p = sin(ϕ) in [-1,1] -> Angle of incidence on boundary reflection
q in [0,length(boundary)] -> Arc length along boundary state is centered at.
ρ -> Density of coherent function, small ρ represents low position uncertainty with
     high momentum uncertainty, large ρ represents the opposite tradeoff. Values around
     1 are usually a good balance.
kn -> k*ni
L -> length of boundary
qp -> q′ the point the state function is evaluated at
mTerms -> Affects number of repetitions of the function along the boundary (similar to m
    quantum number)
"""
function coherentState(p, q, ρ, kn, L, qp; mTerms=5)
    prefactor = (kn/(π*ρ^2))^(1/4)
    sumTerms = [exp(im*p*kn*(qp-q+m*L)-kn/(2ρ^2)*(qp-q+m*L)^2) for m in -mTerms:mTerms]
    prefactor*sum(sumTerms)
end

"""
Calculate the Husimi function from the boundary field data (generally normal derivative)

p = sin(ϕ) in [-1,1] -> Angle of incidence on boundary reflection
q in [0,length(boundary)] -> Arc length along boundary state is centered at.
ρ -> Density of coherent function, small ρ represents low position uncertainty with
     high momentum uncertainty, large ρ represents the opposite tradeoff. Values around
     1 are usually a good balance.
k -> Resonant wavenumber used to calculate field
boundary -> Boundary in question
∂νΨ -> Normal derivative of wavefunction on boundary (NOTE: this is the second part of
       xf as returned from the BoundaryIntegralMethod.residualInverseIteration routine).
       In practice, besides normalization effects the field value (first half of xf) can
       also work well.
qs -> arcLength traversed (q coordinate) at same points as for ∂νΨ
dqs -> local dq/dt at the same points as for ∂νΨ
mTerms -> Affects number of repetitions of the function along the boundary (similar to m
    quantum number)
"""
function husimi(p::Real, q::Real, ρ::Real, k::Real, L::Real, boundary::Boundary, ∂νΨ::Array{Complex{Float64}}, qs::Array{Float64}, dqs::Array{Float64}; mTerms=5)
    @assert length(∂νΨ) == length(qs) == length(dqs)
    kn = k*nInt(boundary)
    sumTerms = [dqs[i]*∂νΨ[i]*coherentState(p, q, ρ, kn, L, qs[i]; mTerms=mTerms) for i in 1:length(∂νΨ)]
    #normalization is not correct, but does it matter?
    1/(2π*abs(kn))*abs(sum(sumTerms))^2
end

"""
Calculate the Husimi distribution for the given boundary, resonant k (kf) and calculated
boundary field (xf), using the momentum resolution specified by pStep.

Returns tuple (qs, ps, hVals)
"""
function husimiDistribution(kf, xf, boundary::Boundary; mTerms=5, pStep=nothing, ρ=1, useField=false)
    ∂νΨ = xf[div(length(xf),2)+1:end]
    if (useField) # Force use of boundary field instead of normal derivative
        ∂νΨ = xf[1:div(length(xf),2)]
    end
    steps = length(∂νΨ)
    stepRange = 0:steps-1
    dqs = arcLen.(Ref(boundary), 2π*stepRange/steps, steps)
    qs = cumsum(vcat([0.0], dqs))
    L = qs[end]
    qs = qs[1:end-1]
    if pStep == nothing
        pStep = 1/steps
    end
    ps = -1:pStep:1
    hVals = @DArray [husimi(p, q, ρ, real(kf), L, boundary, ∂νΨ, qs, dqs; mTerms=mTerms) for p in ps, q in qs]
    (qs, ps, BoundaryIntegralMethod.convertAndClose(hVals))
end
end
