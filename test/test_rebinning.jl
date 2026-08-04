using Test
using Stingray
using Statistics

@testset "Rebinning Utilities" begin

    @testset "rebin_data linear - sum mode snapshot" begin
        # 10 bins with known ascending values, rebin by factor 2
        x = collect(0.0:1.0:9.0)
        y = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        xbin, ybin, _, step = rebin_data(x, y, 2.0; method=:sum, dx=1.0)
        @test length(ybin) == 5
        @test xbin ≈ [1.0, 3.0, 5.0, 7.0, 9.0]
        @test ybin ≈ [3.0, 7.0, 11.0, 15.0, 19.0]
        @test all(step .≈ 2.0)
    end

    @testset "rebin_data linear - mean mode snapshot" begin
        # Same 10 ascending bins, mean should be the average of each pair
        x = collect(0.0:1.0:9.0)
        y = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        xbin, ybin, _, _ = rebin_data(x, y, 2.0; method=:mean, dx=1.0)
        @test xbin ≈ [1.0, 3.0, 5.0, 7.0, 9.0]
        @test ybin ≈ [1.5, 3.5, 5.5, 7.5, 9.5]

        # Also check factor-5 rebinning: 10 bins → 2 bins
        xbin5, ybin5, _, _ = rebin_data(x, y, 5.0; method=:mean, dx=1.0)
        @test xbin5 ≈ [2.5, 7.5]
        @test ybin5 ≈ [3.0, 8.0]  # mean(1:5)=3, mean(6:10)=8
    end

    @testset "rebin_data linear - error propagation" begin
        x = collect(0.0:1.0:9.0)
        y = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        yerr = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        # Sum mode: error = sqrt(e1^2 + e2^2) per bin pair
        _, _, ybinerr_sum, _ = rebin_data(x, y, 2.0; yerr=yerr, method=:sum, dx=1.0)
        @test ybinerr_sum[1] ≈ sqrt(0.1^2 + 0.2^2)  # ≈ 0.2236
        @test ybinerr_sum[2] ≈ sqrt(0.3^2 + 0.4^2)  # = 0.5
        @test ybinerr_sum[3] ≈ sqrt(0.5^2 + 0.6^2)  # ≈ 0.7810

        # Mean mode: error = sqrt(e1^2 + e2^2) / nsamples
        _, _, ybinerr_mean, _ = rebin_data(x, y, 2.0; yerr=yerr, method=:mean, dx=1.0)
        @test ybinerr_mean[1] ≈ sqrt(0.1^2 + 0.2^2) / 2.0
        @test ybinerr_mean[2] ≈ sqrt(0.3^2 + 0.4^2) / 2.0
    end

    @testset "rebin_data linear - complex values" begin
        x = collect(0.0:1.0:9.0)
        yc = Complex{Float64}[1+0.5im, 2+1im, 3+0.5im, 4+1im, 5+0.5im,
                               6+1im, 7+0.5im, 8+1im, 9+0.5im, 10+1im]
        _, ybin_c, _, _ = rebin_data(x, yc, 2.0; method=:mean, dx=1.0)
        @test ybin_c[1] ≈ 1.5 + 0.75im
        @test ybin_c[2] ≈ 3.5 + 0.75im
        @test ybin_c[3] ≈ 5.5 + 0.75im
    end

    @testset "rebin_data linear - resolution guard" begin
        x = collect(0.0:0.01:0.99)
        y = ones(length(x))
        @test_throws ArgumentError rebin_data(x, y, 0.005; dx=0.01)
    end

    @testset "rebin_data linear - non-uniform binning" begin
        # Without specifying dx, should compute from diff(x)
        x = collect(0.0:1.0:9.0)
        y = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        _, ybin, _, _ = rebin_data(x, y, 2.0; method=:mean)
        @test ybin ≈ [1.5, 3.5, 5.5, 7.5, 9.5]
    end

    @testset "rebin_data_log - snapshot" begin
        freq = collect(1.0:1.0:10.0)
        power = Float64[2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        xlog, ylog, _, nsamp = rebin_data_log(freq, power, 0.3; dx=1.0)
        @test length(ylog) < length(power)
        @test all(nsamp .>= 1)
        # Check exact rebinned values
        @test xlog ≈ [1.0, 2.0, 3.5, 5.5, 8.0, 10.0]
        @test ylog ≈ [2.0, 4.0, 7.0, 11.0, 16.0, 20.0]
        @test nsamp == [1, 1, 2, 2, 3, 1]
        # nsamples should be non-decreasing (ignoring edge effects at the last bin)
        @test issorted(nsamp[1:end-1])
    end

    @testset "rebin_data_log - complex values" begin
        freq = collect(1.0:1.0:10.0)
        pc = Complex{Float64}[1+0.5im, 2+1im, 3+0.5im, 4+1im, 5+0.5im,
                               6+1im, 7+0.5im, 8+1im, 9+0.5im, 10+1im]
        xlog, ylog_c, _, _ = rebin_data_log(freq, pc, 0.3; dx=1.0)
        @test eltype(ylog_c) <: Complex
        # First bin is a single sample, should equal the original
        @test ylog_c[1] ≈ 1.0 + 0.5im
        # Third bin averages bins 3 and 4: mean(3+0.5im, 4+1im) = 3.5+0.75im
        @test ylog_c[3] ≈ 3.5 + 0.75im
    end

    @testset "rebin_data_log - input validation" begin
        @test_throws ArgumentError rebin_data_log([1.0, 2.0], [1.0], 0.01)
    end
end
