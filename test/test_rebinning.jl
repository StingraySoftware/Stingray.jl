using Test
using Stingray
using Statistics

@testset "Rebinning Tests" begin
    # Helper to create dummy LightCurve
    function make_lc(times, counts, dt)
        meta = LightCurveMetadata("test", "test", "test", 0.0, (times[1], times[end]), dt, [], Dict{String,Any}())
        exposure = fill(dt, length(times))
        LightCurve{Float64}(times, dt, counts, nothing, exposure, EventProperty{Float64}[], meta, :poisson)
    end

    @testset "LightCurve Rebinning Overload" begin
        dt = 1.0
        times = collect(0:dt:100.0-dt)
        counts = ones(Int, length(times))
        lc = make_lc(times, counts, dt)
        
        # Test rebin(lc, factor)
        factor = 2
        lc_binned = rebin(lc, factor)
        
        @test lc_binned.dt == dt * factor
        @test length(lc_binned.time) == length(lc.time) ÷ factor
        @test sum(lc_binned.counts) == sum(lc.counts) # Invariant: Total counts preserved
        @test lc_binned.metadata.bin_size == dt * factor
    end

    @testset "PowerSpectrum Rebinning (Linear)" begin
        dt = 0.1
        times = collect(0:dt:100.0-dt)
        n_bins = length(times)
        counts = rand(0:10, n_bins) 
        lc = make_lc(times, counts, dt)
        ps = Powerspectrum(lc, norm="leahy")
        
        factor = 4
        ps_binned = rebin(ps, factor)
        
        @test length(ps_binned.freq) == floor(Int, length(ps.freq) / factor)
        @test ps_binned.df ≈ ps.df * factor
        
        # Check monotonicity
        @test all(diff(ps_binned.freq) .> 0)
    end

    @testset "PowerSpectrum Rebinning (Log)" begin
        dt = 0.01
        times = collect(0:dt:100.0-dt)
        counts = rand(0:10, length(times))
        lc = make_lc(times, counts, dt)
        ps = Powerspectrum(lc)
        
        f_res = 0.1
        ps_log = logrebin(ps, f=f_res)
        
        # Check that frequencies are roughly increasing geometrically
        ratios = ps_log.freq[2:end] ./ ps_log.freq[1:end-1]
        # This is a loose check because bins are discrete
        @test all(ratios .>= 1.0)
        
        # Check that we have fewer bins than original
        @test length(ps_log.freq) < length(ps.freq)
        
        # Check monotonicity
        @test all(diff(ps_log.freq) .> 0)
        
        # Check that df is set (not 0.0)
        @test ps_log.df > 0.0
    end
end
