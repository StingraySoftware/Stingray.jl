using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging

TEST_DATA_PATH = joinpath(@__DIR__, "_data")

TEST_DATA_REMOTE = "https://www.star.bristol.ac.uk/fergus/stingray/testdata.tar.gz"
TEST_DATA_SHA256 = "cd5311e5ea1aaf7eef8253a369bd8ba6fbc10c922a22281adc14be61dc6e468c"

# assert the test data is available
if !isdir(TEST_DATA_PATH)
    using Downloads, Tar, CodecZlib, SHA
    dest = TEST_DATA_PATH * ".tar.gz"
    Downloads.download(TEST_DATA_REMOTE, dest)

    open(dest) do tar_gz
        @assert bytes2hex(sha256(tar_gz)) == TEST_DATA_SHA256
        seek(tar_gz, 0)
        tar = GzipDecompressorStream(tar_gz)
        Tar.extract(tar, TEST_DATA_PATH)
    end

    rm(dest)
else
    @info "Test data is present"
end

@testset "GTI" begin
    include("test_gti.jl")
end

@testset "Fourier" begin
    include("test_fourier.jl")
end

@testset "Eventlist" begin
    include("test_events.jl")
end
