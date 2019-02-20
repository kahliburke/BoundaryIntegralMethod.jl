module Utilities
using BoundaryIntegralMethod
using BoundaryIntegralMethod.Husimi
using JLD

export writeWithHusimi

function writeWithHusimi(filename, b, qs, ps, hVals, kf, xf, field)
    boundaryField = [boundaryPoint(b, t) for t in 0.:pi/length(xf):2pi][1:end-1]
    ts = (0:2*pi/length(qs):2pi)[1:end-1]
    jldopen(filename, "w") do file
        write(file, "qs", qs)
        write(file, "ps", convert(Array, ps))
        write(file, "hVals", hVals)
        write(file, "ts", convert(Array,ts))
        write(file, "kf/re", real(kf))
        write(file, "kf/im", imag(kf))
        write(file, "kf/full", kf)
        write(file, "xf/re", real.(xf))
        write(file, "xf/im", imag.(xf))
        write(file, "xf/full", xf)
        write(file, "field/re", real.(field))
        write(file, "field/im", imag.(field))
        write(file, "field/full", field)
        write(file, "boundaryPoints/full", boundaryField)
        write(file, "boundaryPoints/x", [p.x for p in boundaryField])
        write(file, "boundaryPoints/y", [p.y for p in boundaryField])
        write(file, "boundary", b)
    end
end
end
