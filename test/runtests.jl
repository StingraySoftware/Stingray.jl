using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,Plots
include("test_fourier.jl")
include("test_gti.jl")
include("test_events.jl")
include("test_powerspectrum.jl")
include("test_lightcurve.jl")
include("test_recipes.jl")