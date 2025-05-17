using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging
include("test_fourier.jl")
include("test_gti.jl")
include("test_events.jl")
include("test_lightcurve.jl")