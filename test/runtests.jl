using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
using Random
using CairoMakie
using RecipesBase

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

include("test_power_colors.jl")

include("test_rebinning.jl")
include("test_crossspectrum.jl")
include("test_powerspectrum.jl")
include("test_lombscargle.jl")