using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
include("test_fourier.jl")
@testset "GTI" begin
    include("test_gti.jl")
end
@testset "Eventlist" begin
    include("test_events.jl")
end

@testset "lightcurve" begin
    include("test_lightcurve.jl")
end