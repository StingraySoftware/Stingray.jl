using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
using Random

include("test_fourier.jl")
@testset "GTI" begin
    include("test_gti.jl")
end
@testset "Eventlist" begin
    include("test_events.jl")
end

@testset "Synthetic Events Tests" begin
    include("test_missionSupport.jl")
end