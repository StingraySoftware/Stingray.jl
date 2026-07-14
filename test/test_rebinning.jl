using Test
using Stingray
using Statistics

@testset "Rebinning Utilities" begin

    @testset "rebin_data linear - sum mode" begin
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        xbin, ybin, ybinerr, step = rebin_data(x, y, 0.04; method=:sum, dx=0.01)
        @test length(ybin) == 25
        @test all(isapprox.(ybin, 4.0, atol=0.1))
        @test all(step .≈ 4.0)
    end

    @testset "rebin_data linear - mean mode" begin
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        _, ybin, _, _ = rebin_data(x, y, 0.04; method=:mean, dx=0.01)
        @test all(isapprox.(ybin, 1.0, atol=0.01))
    end

    @testset "rebin_data linear - error propagation" begin
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        yerr = fill(2.0, length(x))
        _, _, ybinerr_sum, _ = rebin_data(x, y, 0.04; yerr=yerr, method=:sum, dx=0.01)
        # Sum of 4 errors of 2.0: sqrt(4 * 4) = 4.0
        @test all(isapprox.(ybinerr_sum, 4.0, atol=0.5))

        _, _, ybinerr_mean, _ = rebin_data(x, y, 0.04; yerr=yerr, method=:mean, dx=0.01)
        # Mean: sqrt(4 * 4) / 4 = 1.0
        @test all(isapprox.(ybinerr_mean, 1.0, atol=0.15))
    end

    @testset "rebin_data linear - complex values" begin
        x = collect(0.0:0.01:0.99)
        yc = ones(Complex{Float64}, length(x)) .+ 0.5im
        _, ybin_c, _, _ = rebin_data(x, yc, 0.04; method=:mean, dx=0.01)
        @test all(isapprox.(real.(ybin_c), 1.0, atol=0.01))
        @test all(isapprox.(imag.(ybin_c), 0.5, atol=0.01))
    end

    @testset "rebin_data linear - resolution guard" begin
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        @test_throws ArgumentError rebin_data(x, y, 0.005; dx=0.01)
    end

    @testset "rebin_data linear - non-uniform binning" begin
        # Without specifying dx, should compute from diff(x)
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        _, ybin, _, _ = rebin_data(x, y, 0.04; method=:mean)
        @test all(isapprox.(ybin, 1.0, atol=0.01))
    end

    @testset "rebin_data_log - basic" begin
        freq = collect(1.0:1.0:100.0)
        power = ones(length(freq))
        xlog, ylog, _, nsamp = rebin_data_log(freq, power, 0.01; dx=1.0)
        @test length(ylog) < length(power)
        @test all(nsamp .>= 1)
        @test all(isapprox.(ylog, 1.0, atol=0.01))
    end

    @testset "rebin_data_log - nsamples growth" begin
        freq = collect(1.0:1.0:100.0)
        power = ones(length(freq))
        _, _, _, nsamp = rebin_data_log(freq, power, 0.01; dx=1.0)
        # Higher frequency bins should have more samples
        @test nsamp[end] >= nsamp[1]
    end

    @testset "rebin_data_log - complex values" begin
        freq = collect(1.0:1.0:100.0)
        pc = ones(Complex{Float64}, length(freq)) .+ 0.5im
        _, ylog_c, _, _ = rebin_data_log(freq, pc, 0.01; dx=1.0)
        @test eltype(ylog_c) <: Complex
        @test all(isapprox.(real.(ylog_c), 1.0, atol=0.01))
        @test all(isapprox.(imag.(ylog_c), 0.5, atol=0.01))
    end

    @testset "rebin_data_log - input validation" begin
        @test_throws ArgumentError rebin_data_log([1.0, 2.0], [1.0], 0.01)
    end
end
