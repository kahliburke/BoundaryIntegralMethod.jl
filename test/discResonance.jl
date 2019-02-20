using Test
using BoundaryIntegralMethod
using BoundaryIntegralMethod.DirectSweep


# Create a circular boundary
# Internal n
ni=2.0
# External n
ne=1.0
# Radius
r = 1.0
b = BoundaryIntegralMethod.Circle(r, ni, ne)

kTrue = 13.708675947990717 - 0.0006855557536092232im
k0= 13.7 -.0001im
n=100
(ks, xs) = directSweep(b, n, k0)
convergence = 1e-6
i = 1
(kf, xf) = residualInverseIteration(25, convergence, b, n, ks[i], xs[i,:])
@test abs(kTrue-kf) ≈ 0 atol=10e-10
