using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
using Random

@testset "Fourier" begin
    include("test_fourier/test_fourier.jl")
    include("test_fourier/test_coherence.jl")
    include("test_fourier/test_norm.jl")
end
@testset "GTI" begin
    include("test_gti.jl")
end
@testset "Eventlist" begin
    include("test_events.jl")
end

@testset "lightcurve" begin
    include("test_lightcurve.jl")
end