module DirectSweep
using BoundaryIntegralMethod
using PyCall
using LinearAlgebra

const linalg = PyNULL()
function __init__()
    copy!(linalg, pyimport_conda("scipy.linalg", "scipy"))
end

export directSweep

function directSweep(b::Boundary, n::Int, k0::Complex)
  @info "Doing initial precalculation for k=$k0"
  bpc0 = BoundaryIntegralMethod.doPrecalculation(b, n, k0)
  @info "Calculating A0"
  A0 = A(b, n, k0, bpc0)
  @info "Calculating A'"
  Apr = A′(b, n, k0, bpc0)
  try
    @info "Solving generalized eigenvalue problem"
    (vals, vecs) = linalg.eig(-A0, Apr)
    #TODO: Investigate performance of built in eigensolver
    # aEigen = eigen(-A0, Apr)
    # vals = aEigen.values
    # vecs = aEigen.vectors
    p = sortperm(abs.(vals))
    ks = k0.+vals[p]
    xs = vecs[:,p]'
    filterIdx = findall(k->imag(k) < 0 && imag(k) > -50,ks)
    ks = ks[filterIdx]
    xs = xs[filterIdx,:]
    (ks, xs)
  catch e
    @error "Got error attempting directSweep" exception=e
    rethrow(e)
  end
end
end
