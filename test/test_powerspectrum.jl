using Test
using Stingray

@testset "Powerspectrum" begin

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

end

@testset "DynamicalPowerspectrum" begin

    test_gti = [0.0 100.0]

    function make_eventlist_with_gti(t, gti_matrix)
        meta = FITSMetadata(
            "[test]", 1, nothing,
            Dict{String,Vector}(), Dict{String,Any}(),
            gti_matrix, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing
        )
        EventList(t, nothing, meta)
    end

    # Create test data: 100s of events, ~100 counts/s
    rng = MersenneTwister(20240202)
    event_times = sort(rand(rng, Uniform(0.0, 100.0), 10000))
    ev = make_eventlist_with_gti(event_times, test_gti)
    dt = 0.01
    segment_size = 10.0

    @testset "EventList constructor" begin
        dps = DynamicalPowerspectrum(ev, segment_size, dt; norm="frac", silent=true)

        @test dps isa DynamicalPowerspectrum{Float64}
        @test dps isa AbstractPowerspectrum
        @test dps.type == "dynamical_powerspectrum"
        @test dps.norm == "frac"
        @test length(dps.freq) > 0
        @test length(dps.time) > 0
        @test dps.df > 0
        @test dps.dt == segment_size
        @test dps.m == 1
        @test dps.nphots > 0
        @test dps.meanrate > 0
        @test eltype(dps.dyn_ps) <: Real
    end

    @testset "Matrix shape" begin
        dps = DynamicalPowerspectrum(ev, segment_size, dt; norm="frac", silent=true)

        n_freq, n_time = size(dps.dyn_ps)
        @test n_freq == length(dps.freq)
        @test n_time == length(dps.time)
        # 100s / 10s = 10 segments
        @test n_time == 10
    end

    @testset "Bad args - short segment" begin
        @test_throws ArgumentError DynamicalPowerspectrum(
            ev, 0.001, dt; silent=true)
    end

    @testset "LightCurve constructor" begin
        lc = create_lightcurve(ev, dt)

        dps = DynamicalPowerspectrum(lc, segment_size; norm="frac", silent=true)

        @test dps isa DynamicalPowerspectrum{Float64}
        @test dps.type == "dynamical_powerspectrum"
        @test size(dps.dyn_ps, 2) == length(dps.time)
    end

    @testset "LightCurve bad args - long segment" begin
        lc = create_lightcurve(ev, dt)

        @test_throws ArgumentError DynamicalPowerspectrum(
            lc, 200.0; silent=true)
    end

    @testset "Base.show" begin
        dps = DynamicalPowerspectrum(ev, segment_size, dt; norm="frac", silent=true)
        buf = IOBuffer()
        show(buf, dps)
        s = String(take!(buf))
        @test occursin("DynamicalPowerspectrum", s)
        @test occursin("frac", s)
        @test occursin("freq bins", s)
        @test occursin("time segments", s)
    end

    # Method tests
    dps = DynamicalPowerspectrum(ev, segment_size, dt; norm="frac", silent=true)

    @testset "rebin_frequency" begin
        df_new = dps.df * 5
        dps_rb = rebin_frequency(dps, df_new)
        @test dps_rb isa DynamicalPowerspectrum{Float64}
        @test length(dps_rb.freq) < length(dps.freq)
        @test size(dps_rb.dyn_ps, 2) == size(dps.dyn_ps, 2)  # time unchanged
    end

    @testset "rebin_time" begin
        dt_new = dps.dt * 2
        dps_rb = rebin_time(dps, dt_new)
        @test dps_rb isa DynamicalPowerspectrum{Float64}
        @test length(dps_rb.time) < length(dps.time)
        @test size(dps_rb.dyn_ps, 1) == size(dps.dyn_ps, 1)  # freq unchanged
    end

    @testset "rebin_by_n_intervals" begin
        dps_rb = rebin_by_n_intervals(dps, 2)
        @test dps_rb isa DynamicalPowerspectrum{Float64}
        @test size(dps_rb.dyn_ps, 2) == size(dps.dyn_ps, 2) ÷ 2
        @test dps_rb.m == dps.m * 2
    end

    @testset "trace_maximum" begin
        idx = trace_maximum(dps)
        @test length(idx) == size(dps.dyn_ps, 2)
        @test all(1 .<= idx .<= length(dps.freq))
    end

    @testset "trace_maximum with bounds" begin
        mid_freq = dps.freq[length(dps.freq) ÷ 2]
        idx = trace_maximum(dps; min_freq=mid_freq)
        @test all(dps.freq[idx] .>= mid_freq)
    end

    @testset "shift_and_add" begin
        f0_list = dps.freq[argmax(dps.dyn_ps[:, j]) for j in 1:size(dps.dyn_ps, 2)]
        result = shift_and_add(dps, f0_list; nbins=50)
        @test haskey(result, :freq)
        @test haskey(result, :power)
        @test haskey(result, :m)
        @test length(result.freq) == 50
        @test length(result.power) == 50
    end

    @testset "compute_rms" begin
        fmin = dps.freq[2]
        fmax = dps.freq[end-1]
        result = compute_rms(dps, fmin, fmax)
        @test length(result.rms) == size(dps.dyn_ps, 2)
        @test length(result.rms_err) == size(dps.dyn_ps, 2)
        @test all(result.rms .>= 0)
    end

    @testset "Normalizations" begin
        for norm in ["frac", "abs", "leahy", "none"]
            dps_n = DynamicalPowerspectrum(ev, segment_size, dt; norm=norm, silent=true)
            @test dps_n.norm == norm
            @test size(dps_n.dyn_ps, 2) > 0
        end
    end

end
