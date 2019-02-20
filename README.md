# BoundaryIntegralMethod

This package implements a boundary integral method approach for solving the Helmholz equation
and finding scattering resonances in dielectric microcavities.

This package includes a set of command line utilities that can be applied to implement
the common steps involved in finding resonances and calculating the electric field in
a region inside and around the microcavity.

The package can also be utilized within a user written self contained script or
within a Jupyter/IJulia notebook.

Several example boundary shapes are implemented and it is easy to add your own boundary
parameterization.

This implementation is based on the derivation by P. Heider<sup id="a1">[1](#f1)</sup> and is also utilized
in an upcoming publication in preparation on folded-chaotic whispering gallery modes.

# Step 1

Use direct sweep algorithm to find candidate resonances close to a starting k

    runDirectSweep.jl --simFile=circle_simple.jl --n=100 --k0="8.9101-0.264im"

The output will be written to an 'output' directory beneath the working directory. This will contain
a file that contains the candidate k values and fields on the boundary for the candidate resonances
found by the direct method. Mathematica (or another tool) could be used to view the resulting HDF5
file and choose which candidate to use.

Run the runDirectSweep.jl script without any arguments to see the available options.

# Step 2

Use residual inverse iteration algorithm to converge on a specific resonance as determined in step 1

    resonantModeSearch.jl --fieldFile=output/sweep-circle_simple-100-8.910171115740056-0.26429086620242176im.h5 --fieldIdx=1 --simFile=circle_simple.jl

In the command above, we set the fieldIdx to 1 to use the first candidate in the fieldFile. This would
correspond to the candidate closest (in terms of absolute distance in the complex plane) to the starting points
specified in step 1.

The output of this process, if the algorithm converges, will be a file in the output directory corresponding
to the final calculated field values on the boundary (boundaryField-xxx.h5), as well as a total field file calculated using the
Green's function from the boundary field (field-xxx.h5).

This field file contains the real and imaginary parts of the field and could be viewed in Mathematica or
other suitable software.

Run the resonantModeSearch.jl script without any arguments to see the available options.

# Step 3

If you want to produce a different field calculation, with a different grid spacing or range, a script is provided. For example,

    fieldFromBoundary.jl --boundaryFile=output/boundaryField-circle_simple-100-8.910171115740056-0.26429086620242176im.h5 --simFile=circle_simple.jl --dx=.005 --dy=.005

Run the fieldFromBoundary.jl script without any arguments to see the available options.

# Step 4

Visualize output by producing a PNG file with the command:

    plotFieldImage.jl --fieldFile=output/field-circle_simple-100-8.910171115740056-0.26429086620242176im.h5 --clamp=1

The clamp option is used to reduce the dynamic range of the field output, it defaults to 1. When producing
the plot, the max absolute value of the field will be shown to assist in choosing the value.

<b id="f1">1</b> P. Heider, Computers & Mathematics With Applications 60, 1620 (2010), ISSN 0898-1221. [↩](#a1)
