using Test
using Stingray
using FFTW, Distributions, Statistics, StatsBase, Metadata, HDF5

@testset "Stingray.jl Tests" begin
    include("test_fourier.jl")
    include("test_gti.jl")
    include("test_logging.jl")
end
