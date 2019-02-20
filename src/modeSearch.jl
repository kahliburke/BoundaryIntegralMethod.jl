using ProgressMeter
using JLD
using LinearAlgebra

function residualInverseIteration(iMax::Int, ϵ::Float64, b::Boundary, n::Int, k0::Complex{Float64}, x0::Array{Complex{Float64},1}; useWiersigKShift::Bool = false)
  residualInverseIteration(iMax, ϵ, b, n, k0, x0, nothing, 0; useWiersigKShift = useWiersigKShift)
end

function residualInverseIteration(iMax::Int, ϵ::Float64, b::Boundary, n::Int, k0::Complex{Float64}, x0::Array{Complex{Float64},1}, progressMeter, iter; useWiersigKShift::Bool = false)
  if progressMeter == nothing
    progressMeter = ProgressThresh(ϵ, "Finding resonance:")
  end

  @debug "Starting iteration, attempt to converge k with ϵ = $ϵ"

  numChoices = min(n,20)
  # Choose from a random large value in x0, don't use the same index each iteration
  # Not sure if this helps in any way
  xMaxInd = sortperm(-abs.(x0[1:n]))[rand(1:numChoices)]
  # xMaxInd = indmax(abs(x0))

  # Construct normalization vector
  e = [j == xMaxInd ? 1.0+0.0im : 0.0+0.0im for j in 1:4n]
  # Normalize incoming x0
  x0 = x0 ./ (e' * x0)[1]
  @debug "max idx = $(xMaxInd), x = $(x0[xMaxInd]), e=$(Array(e)), x = $(x0)"

  @debug "Doing initial precalculation for k=$k0"
  bpc0 = doPrecalculation(b, n, k0)
  @debug "Calculating A matrix"
  A0 = A(b, n, k0, bpc0)
  @debug "Inverting A matrix"
  A0inv = inv(A0)

  xl = x0
  kl = k0
  ks = [kl]

  for l in 1:iMax
    iter += 1
    @debug "l = $l, doing precalculation for k=$kl"
    bpcl = doPrecalculation(b, n, kl)
    @debug "Calculating A' derivative matrix"
    A′l = A′(b, n, kl, bpcl)
    @debug "Calculating A(kl) matrix"
    Akl = A(b, n, kl, bpcl)

    @debug "e' * (A0inv*Akl*xl)" dote = e' * (A0inv*Akl*xl)
    @debug "e' * (A0inv*A′l*xl)" dotA = e' * (A0inv*A′l*xl)

    if useWiersigKShift
      # From Wiersig's paper eq (37), this seems to work slightly more slowly and less
      # accurately than Heider's version. I'm curious if there are ever cases
      # where it could work better?
      Aklinv = inv(Akl)
      dk = - 1/(trace(Aklinv*A′l))
      dkHeider = - ((e' * (A0inv*Akl*xl)) / (e' * (A0inv*A′l*xl)))[1]
      @debug "Heider's dk would have been = $(dkHeider)"
    else
      dk = - ((e' * (A0inv*Akl*xl)) / (e' * (A0inv*A′l*xl)))[1]
    end
    @debug "dk: $(dk)"
    kNext = kl+dk
    @debug "k(l+1): $kNext"
    push!(ks, kNext)
    @debug "Doing l+1 precalculation for k=$kNext"
    bpcNext = doPrecalculation(b, n, kNext)
    @debug "Calculating A matrix"
    Anext = A(b, n, kNext, bpcNext)
    # Calculate residuum
    rl = Anext*xl
    # Use that to determine dx
    # Had this from Heider's paper but since we've already calcuated A0inv, seems like
    # just as good to use it?
    # dxl = A0 \ rl
    dxl = A0inv*rl
    @debug "|dx_l| = $(norm(dxl))"

    #Adjust x vector
    xNextBar = xl - dxl
    xNext = xNextBar ./ ((e' * xNextBar)[1])

    @debug "Checking |A(kl).x_l+1| ~= 0"
    beforeShiftQuality = norm(Anext*xl)/norm(xl)
    afterShiftQuality = norm(Anext*xNext)/norm(xNext)
    @debug "|A(kl+1).x_l|/|x_l|=$beforeShiftQuality"
    @debug "|A(kl+1).x_l+1|/|x_l+1|=$afterShiftQuality"

    kl = kNext
    xl = xNext
    ProgressMeter.update!(progressMeter, abs.(ks[end]-ks[end-1]), showvalues = [
      (:iter,iter),
      (:dk,dk),
      (:kl,kl),
      ("|A(kl+1).x_l|/|x_l|", beforeShiftQuality),
      ("|A(kl+1).x_l+1|/|x_l+1|", afterShiftQuality)
    ])
    if !correctorCondition(ks)
      break
    end
  end

  if length(ks) > 1 && abs.(ks[end]-ks[end-1]) < ϵ
    @debug "k converged to $kl"
    # Clear memory
    bpc0 = 0
    A0 = 0
    A0inv = 0
    x0 = 0
    bpcF = doPrecalculation(b, n, kl)
    Af = A(b, n, kl, bpcF)
    fEigen = eigen(Af)
    fVals = fEigen.values
    fVecs = fEigen.vectors
    perm = sortperm(abs.(fVals))
    @info "Min eigenvalue = $(fVals[perm[1]])"
    xl = Array(xl)
    xfDiff = norm(fVecs[:,perm[1]]-Array(xl))
    @info "Diff in field = $xfDiff"
    (kl, Array(xl))
  elseif (iter < iMax)
    @debug "Restart RII with shift $kl and current xl"
    try
      # Clear memory
      bpc0 = 0
      A0 = 0
      A0inv = 0
      x0 = 0
      residualInverseIteration(iMax, ϵ, b, n, kl, xl, progressMeter, iter; useWiersigKShift = useWiersigKShift)
    catch e
      @error ("Got exception in RII") exception=e
      error("k did not converge after $iMax iterations, last kl = $kl")
    end
  else
    error("k did not converge after $iMax iterations, last kl = $kl")
  end
end

function correctorCondition(ks)
  if size(ks)[1] == 1
    true
  elseif abs.((ks[end]-ks[end-1])/ks[end]) < 1/10
    @debug "Hit corrector condition, kl=$(ks[end])"
    false
  else
    true
  end
end
