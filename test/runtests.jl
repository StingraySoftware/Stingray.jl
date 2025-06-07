using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
include("test_fourier.jl")
include("test_gti.jl")
@testset "Eventlist" begin
    include("test_events.jl")
end

@testset "LightCurve" begin
    include("test_lightcurve.jl")
end