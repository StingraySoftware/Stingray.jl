using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging ,LinearAlgebra
using CFITSIO
include("test_fourier.jl")
include("test_gti.jl")
@testset verbose=true "Eventlist" begin
    include("test_events.jl")
end
@testset verbose=true "Synthetic Events Tests" begin
    include("test_missionSupport.jl")
end