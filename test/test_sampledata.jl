using Test
using DelimitedFiles
using Downloads

# Define the URL of the file later on when we we want to import datsets from stingray software we can do it in here and remove this function this just for following the guidelines of importing
const LC_FILE_URL = "https://raw.githubusercontent.com/StingraySoftware/stingray/main/stingray/datasets/lc_sample.txt"
const LC_FILE_PATH = joinpath(@__DIR__, "lc_sample.txt")
# here i have created function that downloads the datasets from the provided url for testing purposes

if !isfile(LC_FILE_PATH)
    Downloads.download(LC_FILE_URL, LC_FILE_PATH)
end

# Including sampledata.jl
include("../src/sampledata.jl")

# Constants
const FIRST_LINE = (1.109110400703125000e08, 4.120000061392784119e03)
const LAST_LINE = (1.109138078203125000e08, 2.619200039029121399e03)
const FILE_LENGTH = 22143

@testset "Test Sample Data" begin
    @testset "Test File Exists" begin
        @test isfile(LC_FILE_PATH) == true
        lc = sample_data()
        @test !isempty(lc.time)
        @test !isempty(lc.counts)
    end

    @testset "Test File First Line" begin
        lc = sample_data()
        @test lc.time[1] ≈ FIRST_LINE[1]
        @test lc.counts[1] ≈ FIRST_LINE[2]
    end

    @testset "Test File Last Line" begin
        lc = sample_data()
        @test !isempty(lc.time)
        @assert !isempty(lc.time) "LightCurve data is empty"
        @test lc.time[end] ≈ LAST_LINE[1]
        @test lc.counts[end] ≈ LAST_LINE[2]
    end

    @testset "Test File Length" begin
        lc = sample_data()
        @test length(lc.time) == FILE_LENGTH
        @test length(lc.counts) == FILE_LENGTH
    end
end
