using Test
using Stingray

@testset "Powerspectrum" begin

    function make_eventlist_with_gti(t, gti_matrix)
        meta = FITSMetadata(
            "[test]", 1, nothing,
            Dict{String,Vector}(), Dict{String,Any}(),
            gti_matrix, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing
        )
        EventList(t, nothing, meta)
    end

    @testset "Direct struct initialization" begin

        freq = [1.0, 2.0, 3.0]
        power = [1.5, 2.0, 1.8]
        power_err = [0.15, 0.2, 0.18]
        unnorm = [150.0, 200.0, 180.0]
        unnorm_err = [15.0, 20.0, 18.0]
        test_gti = [0.0 10.0]

        ps = Powerspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            1.0, 0.0001, 100000, 1, 1,  # df, dt, n, m, k
            1000.0,  # nphots
            "leahy",  # norm
            test_gti, 10.0,  # gti, segment_size
            nothing,  # variance
            "poisson", "powerspectrum"  # err_dist, type
        )

        @test ps isa Powerspectrum{Float64}
        @test ps isa AbstractPowerspectrum
        @test ps isa AbstractCrossspectrum  # subtype check
        @test length(ps.freq) == 3
        @test ps.norm == "leahy"
        @test ps.type == "powerspectrum"
        @test ps.m == 1
        @test ps.k == 1
        @test isnothing(ps.variance)
        @test ps.err_dist == "poisson"
        @test eltype(ps.power) <: Real
    end

    @testset "Base.show" begin
        freq = [1.0, 2.0, 3.0]
        power = [1.5, 2.0, 1.8]
        test_gti = [0.0 10.0]
        ps = Powerspectrum{Float64}(
            freq, power, power, power, power,
            1.0, 0.0001, 100000, 1, 1,
            1000.0, "leahy", test_gti, 10.0, nothing,
            "poisson", "powerspectrum"
        )
        buf = IOBuffer()
        show(buf, ps)
        s = String(take!(buf))
        @test occursin("Powerspectrum", s)
        @test occursin("leahy", s)
        @test occursin("1 segments", s)
    end

    @testset "EventList constructor" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        rng = MersenneTwister(20150907)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        ev = EventList(times1, nothing, meta)

        ps = Powerspectrum(ev, 10.0, 0.001; norm="frac", silent=true)

        @test ps isa Powerspectrum{Float64}
        @test ps.norm == "frac"
        @test ps.type == "powerspectrum"
        @test length(ps.freq) > 0
        @test ps.df > 0
        @test ps.m >= 1
        @test ps.k == 1
        @test ps.nphots > 0
        @test ps.err_dist == "poisson"
        @test eltype(ps.power) <: Real
    end

    @testset "EventList constructor normalizations" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        rng = MersenneTwister(42)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        ev = EventList(times1, nothing, meta)

        for norm in ["frac", "abs", "leahy", "none"]
            ps = Powerspectrum(ev, 10.0, 0.001; norm=norm, silent=true)
            @test ps.norm == norm
            @test length(ps.freq) > 0
        end
    end

    @testset "power_err calculation" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        rng = MersenneTwister(123)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        ev = EventList(times1, nothing, meta)

        ps = Powerspectrum(ev, 10.0, 0.001; norm="frac", silent=true)
        expected_err = ps.power ./ sqrt(ps.m)
        @test ps.power_err ≈ expected_err
    end

    @testset "LightCurve constructor" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        rng = MersenneTwister(20150907)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        ev = EventList(times1, nothing, meta)

        lc = create_lightcurve(ev, 0.001)

        ps = Powerspectrum(lc, 5.0; norm="frac", silent=true)

        @test ps isa Powerspectrum{Float64}
        @test ps.norm == "frac"
        @test ps.type == "powerspectrum"
        @test length(ps.freq) > 0
        @test ps.m >= 1
        @test ps.k == 1
        @test ps.err_dist == "poisson"
    end

    @testset "Leahy normalization Poisson noise" begin
        # Simulate Poisson noise
        rng = MersenneTwister(42)
        rate = 100.0  # counts/s
        duration = 100.0
        dt = 0.01
        times = sort(rand(rng, Uniform(0.0, duration), round(Int, rate * duration)))
        
        test_gti = [0.0 duration]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        ev = EventList(times, nothing, meta)

        ps = Powerspectrum(ev, duration, dt; norm="leahy", silent=true)
        
        # In Leahy norm, Poisson noise should have mean power ≈ 2.0
        # Exclude the DC bin (first bin) if necessary, but here we can check the whole array
        @test isapprox(mean(ps.power[2:end]), 2.0, atol=0.1)
    end

    @testset "Fractional RMS and abs normalization" begin
        rng = MersenneTwister(42)
        rate = 100.0  # counts/s
        duration = 100.0
        dt = 0.01
        times = sort(rand(rng, Uniform(0.0, duration), round(Int, rate * duration)))
        times[1] = 0.0
        times[end] = duration
        
        test_gti = [0.0 duration]
        meta = FITSMetadata(
            "[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(),
            test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        )
        ev = EventList(times, nothing, meta)
        lc = create_lightcurve(ev, dt)

        # test fractional RMS normalization
        ps_frac = Powerspectrum(lc, duration; norm="frac", silent=true)
        # fractional RMS squared should approximately equal var(counts)/mean(counts)^2
        rms_sq = sum(ps_frac.power[2:end] .* ps_frac.df)
        expected_rms_sq = var(lc.counts) / mean(lc.counts)^2
        # Poisson noise has a specific level, we can check if they are of the same order
        # or exactly equal depending on the implementation details.
        @test isapprox(rms_sq, expected_rms_sq, rtol=0.2)

        # test abs normalization
        ps_abs = Powerspectrum(lc, duration; norm="abs", silent=true)
        # For abs norm, Poisson noise level should be 2 * countrate
        countrate = ps_abs.nphots / ps_abs.segment_size
        @test isapprox(mean(ps_abs.power[2:end]), 2.0 * countrate, atol=0.1 * countrate)
    end

    @testset "rebin - linear snapshot" begin
        # Build a Powerspectrum with known values (10 bins) for exact checks
        freq = collect(1.0:1.0:10.0)
        power = Float64[2.1, 1.8, 2.3, 1.9, 2.0, 2.2, 1.7, 2.4, 2.0, 1.6]
        power_err = power ./ sqrt(1)
        unnorm = power .* 100.0
        unnorm_err = power_err .* 100.0
        test_gti = [0.0 10.0]

        ps = Powerspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            1.0, 0.001, 10000, 1, 1,
            1000.0, "leahy", test_gti, 10.0, nothing, "poisson", "powerspectrum"
        )

        # Rebin by factor 2 (10 → 5 bins)
        ps_rb = rebin(ps, 2.0; method=:mean)
        @test ps_rb.freq ≈ [2.0, 4.0, 6.0, 8.0, 10.0]
        @test ps_rb.power ≈ [1.95, 2.1, 2.1, 2.05, 1.8]  atol=1e-10
        @test ps_rb.df ≈ 2.0
        @test ps_rb.k == 2
        @test ps_rb.m == 2
        @test eltype(ps_rb.power) <: Real

        # Rebin by factor 5 (10 → 2 bins)
        ps_rb5 = rebin(ps, 5.0; method=:mean)
        @test ps_rb5.freq ≈ [3.5, 8.5]  atol=1e-10
        @test ps_rb5.power ≈ [2.02, 1.98]  atol=1e-10
        @test ps_rb5.k == 5
    end

    @testset "rebin_log - snapshot" begin
        freq = collect(1.0:1.0:10.0)
        power = Float64[2.1, 1.8, 2.3, 1.9, 2.0, 2.2, 1.7, 2.4, 2.0, 1.6]
        power_err = power ./ sqrt(1)
        unnorm = power .* 100.0
        unnorm_err = power_err .* 100.0
        test_gti = [0.0 10.0]

        ps = Powerspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            1.0, 0.001, 10000, 1, 1,
            1000.0, "leahy", test_gti, 10.0, nothing, "poisson", "powerspectrum"
        )

        ps_log = rebin_log(ps, 0.3)
        @test ps_log.freq ≈ [1.0, 2.0, 3.5, 5.5, 8.0, 10.0]
        # 5th bin averages power at freq 7,8,9 → values 1.7, 2.4, 2.0
        @test ps_log.power ≈ [2.1, 1.8, 2.1, 2.1, (1.7+2.4+2.0)/3, 1.6]  atol=1e-10
        @test ps_log.k == [1, 1, 2, 2, 3, 1]
        @test ps_log.m == [1, 1, 2, 2, 3, 1]
        @test eltype(ps_log.power) <: Real
    end

    @testset "compute_rms - Poisson noise" begin
        rng = MersenneTwister(42)
        rate = 100.0  # counts/s
        duration = 100.0
        dt = 0.01
        times = sort(rand(rng, Uniform(0.0, duration), round(Int, rate * duration)))

        test_gti = [0.0 duration]
        ev_rms = make_eventlist_with_gti(times, test_gti)
        ps_rms = Powerspectrum(ev_rms, duration, dt; norm="frac", silent=true)

        rms, rms_err = compute_rms(ps_rms, ps_rms.freq[2], ps_rms.freq[end])

        # For pure Poisson noise, the RMS should be small and positive
        @test rms > 0
        @test rms_err > 0
        # Check RMS is finite and of a physically reasonable magnitude
        # (for Poisson noise at 100 ct/s, fractional RMS should be small)
        @test rms < 1.0
        @test isapprox(rms, 0.077, atol=0.05)
    end

end

