using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, Metadata, HDF5

include("../src/conftest.jl") 

@info "Running tests for Stingray.jl"

@testset "Stingray.jl Tests" begin
    include("test_fourier.jl")
    include("test_gti.jl")
end
