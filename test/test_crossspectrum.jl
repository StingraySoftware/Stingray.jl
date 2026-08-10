using Test
using Stingray

@testset "Crossspectrum" begin

    test_gti = [0.0 10.0]

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
        power = Complex{Float64}.([1.0+0.5im, 2.0+1.0im, 3.0+0.2im])
        power_err = Complex{Float64}.([0.1+0.05im, 0.2+0.1im, 0.3+0.02im])
        unnorm = Complex{Float64}.([10.0+5.0im, 20.0+10.0im, 30.0+2.0im])
        unnorm_err = Complex{Float64}.([1.0+0.5im, 2.0+1.0im, 3.0+0.2im])
        test_gti = [0.0 10.0]

        cs = Crossspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            nothing, nothing,  # pds1, pds2
            1.0, 0.001, 10000, 1, 1,  # df, dt, n, m, k
            100.0, 120.0, 80.0,  # nphots, nphots1, nphots2
            12.0, 8.0,  # countrate1, countrate2
            "frac", "all", false,  # norm, power_type, fullspec
            test_gti, 10.0,  # gti, segment_size
            "poisson", "crossspectrum"  # err_dist, type
        )

        @test cs isa Crossspectrum{Float64}
        @test cs isa AbstractCrossspectrum
        @test length(cs.freq) == 3
        @test cs.norm == "frac"
        @test cs.type == "crossspectrum"
        @test cs.m == 1
        @test cs.k == 1
        @test cs.fullspec == false
        @test cs.err_dist == "poisson"
        @test eltype(cs.power) <: Complex
    end

    @testset "Base.show" begin
        freq = [1.0, 2.0, 3.0]
        power = Complex{Float64}.([1.0, 2.0, 3.0])
        cs = Crossspectrum{Float64}(
            freq, power, power, power, power,
            nothing, nothing,
            1.0, 0.001, 10000, 1, 1,
            100.0, 120.0, 80.0,
            12.0, 8.0,
            "leahy", "all", false,
            test_gti, 10.0,
            "poisson", "crossspectrum"
        )
        buf = IOBuffer()
        show(buf, cs)
        s = String(take!(buf))
        @test occursin("Crossspectrum", s)
        @test occursin("leahy", s)
        @test occursin("1 segments", s)
    end

    @testset "EventList constructor" begin
        rng = MersenneTwister(20150907)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = make_eventlist_with_gti(times1, test_gti)
        ev2 = make_eventlist_with_gti(times2, test_gti)

        cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm="frac", silent=true)

        @test cs isa Crossspectrum{Float64}
        @test cs.norm == "frac"
        @test cs.type == "crossspectrum"
        @test length(cs.freq) > 0
        @test cs.df > 0
        @test cs.m >= 1
        @test cs.k == 1
        @test cs.nphots1 > 0
        @test cs.nphots2 > 0
        @test cs.err_dist == "poisson"
        @test eltype(cs.power) <: Complex
    end

    @testset "EventList constructor normalizations" begin
        rng = MersenneTwister(42)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = make_eventlist_with_gti(times1, test_gti)
        ev2 = make_eventlist_with_gti(times2, test_gti)

        for norm in ["frac", "abs", "leahy", "none"]
            cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm=norm, silent=true)
            @test cs.norm == norm
            @test length(cs.freq) > 0
        end
    end

    @testset "Complex power type" begin
        rng = MersenneTwister(77)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = make_eventlist_with_gti(times1, test_gti)
        ev2 = make_eventlist_with_gti(times2, test_gti)

        cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm="frac", silent=true)
        @test eltype(cs.power) <: Complex
        @test eltype(cs.unnorm_power) <: Complex
        @test eltype(cs.power_err) <: Complex
    end

    @testset "power_err calculation" begin
        rng = MersenneTwister(123)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = make_eventlist_with_gti(times1, test_gti)
        ev2 = make_eventlist_with_gti(times2, test_gti)

        cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm="frac", silent=true)
        expected_err = cs.power ./ sqrt(cs.m)
        @test cs.power_err ≈ expected_err
    end

    @testset "LightCurve constructor" begin
        rng = MersenneTwister(20150907)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = make_eventlist_with_gti(times1, test_gti)
        ev2 = make_eventlist_with_gti(times2, test_gti)

        lc1 = create_lightcurve(ev1, 0.001)
        lc2 = create_lightcurve(ev2, 0.001)

        cs = Crossspectrum(lc1, lc2, 5.0; norm="frac", silent=true)

        @test cs isa Crossspectrum{Float64}
        @test cs.norm == "frac"
        @test cs.type == "crossspectrum"
        @test length(cs.freq) > 0
        @test cs.m >= 1
        @test cs.k == 1
    end

end

@testset "DynamicalCrossspectrum" begin

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
    rng = MersenneTwister(20240101)
    times1 = sort(rand(rng, Uniform(0.0, 100.0), 10000))
    times2 = sort(rand(rng, Uniform(0.0, 100.0), 8000))
    ev1 = make_eventlist_with_gti(times1, test_gti)
    ev2 = make_eventlist_with_gti(times2, test_gti)
    dt = 0.01
    segment_size = 10.0

    @testset "EventList constructor" begin
        dcs = DynamicalCrossspectrum(ev1, ev2, segment_size, dt; norm="frac", silent=true)

        @test dcs isa DynamicalCrossspectrum{Float64}
        @test dcs isa AbstractCrossspectrum
        @test dcs.type == "dynamical_crossspectrum"
        @test dcs.norm == "frac"
        @test length(dcs.freq) > 0
        @test length(dcs.time) > 0
        @test dcs.df > 0
        @test dcs.dt == segment_size
        @test dcs.m == 1
        @test dcs.nphots1 > 0
        @test dcs.nphots2 > 0
        @test eltype(dcs.dyn_ps) <: Complex
    end

    @testset "Matrix shape" begin
        dcs = DynamicalCrossspectrum(ev1, ev2, segment_size, dt; norm="frac", silent=true)

        n_freq, n_time = size(dcs.dyn_ps)
        @test n_freq == length(dcs.freq)
        @test n_time == length(dcs.time)
        # 100s / 10s = 10 segments
        @test n_time == 10
    end

    @testset "Bad args - short segment" begin
        @test_throws ArgumentError DynamicalCrossspectrum(
            ev1, ev2, 0.001, dt; silent=true)
    end

    @testset "LightCurve constructor" begin
        lc1 = create_lightcurve(ev1, dt)
        lc2 = create_lightcurve(ev2, dt)

        dcs = DynamicalCrossspectrum(lc1, lc2, segment_size; norm="frac", silent=true)

        @test dcs isa DynamicalCrossspectrum{Float64}
        @test dcs.type == "dynamical_crossspectrum"
        @test size(dcs.dyn_ps, 2) == length(dcs.time)
    end

    @testset "LightCurve bad args - long segment" begin
        lc1 = create_lightcurve(ev1, dt)
        lc2 = create_lightcurve(ev2, dt)

        @test_throws ArgumentError DynamicalCrossspectrum(
            lc1, lc2, 200.0; silent=true)
    end

    @testset "Base.show" begin
        dcs = DynamicalCrossspectrum(ev1, ev2, segment_size, dt; norm="frac", silent=true)
        buf = IOBuffer()
        show(buf, dcs)
        s = String(take!(buf))
        @test occursin("DynamicalCrossspectrum", s)
        @test occursin("frac", s)
        @test occursin("freq bins", s)
        @test occursin("time segments", s)
    end

    # Method tests
    dcs = DynamicalCrossspectrum(ev1, ev2, segment_size, dt; norm="frac", silent=true)

    @testset "rebin_frequency" begin
        df_new = dcs.df * 5
        dcs_rb = rebin_frequency(dcs, df_new)
        @test dcs_rb isa DynamicalCrossspectrum{Float64}
        @test length(dcs_rb.freq) < length(dcs.freq)
        @test size(dcs_rb.dyn_ps, 2) == size(dcs.dyn_ps, 2)  # time unchanged
    end

    @testset "rebin_time" begin
        dt_new = dcs.dt * 2
        dcs_rb = rebin_time(dcs, dt_new)
        @test dcs_rb isa DynamicalCrossspectrum{Float64}
        @test length(dcs_rb.time) < length(dcs.time)
        @test size(dcs_rb.dyn_ps, 1) == size(dcs.dyn_ps, 1)  # freq unchanged
    end

    @testset "rebin_by_n_intervals" begin
        dcs_rb = rebin_by_n_intervals(dcs, 2)
        @test dcs_rb isa DynamicalCrossspectrum{Float64}
        @test size(dcs_rb.dyn_ps, 2) == size(dcs.dyn_ps, 2) ÷ 2
        @test dcs_rb.m == dcs.m * 2
    end

    @testset "trace_maximum" begin
        idx = trace_maximum(dcs)
        @test length(idx) == size(dcs.dyn_ps, 2)
        @test all(1 .<= idx .<= length(dcs.freq))
    end

    @testset "trace_maximum with bounds" begin
        mid_freq = dcs.freq[length(dcs.freq) ÷ 2]
        idx = trace_maximum(dcs; min_freq=mid_freq)
        @test all(dcs.freq[idx] .>= mid_freq)
    end

    @testset "shift_and_add" begin
        f0_list = dcs.freq[argmax(abs.(dcs.dyn_ps[:, j])) for j in 1:size(dcs.dyn_ps, 2)]
        result = shift_and_add(dcs, f0_list; nbins=50)
        @test haskey(result, :freq)
        @test haskey(result, :power)
        @test haskey(result, :m)
        @test length(result.freq) == 50
        @test length(result.power) == 50
    end

    @testset "compute_rms" begin
        fmin = dcs.freq[2]
        fmax = dcs.freq[end-1]
        result = compute_rms(dcs, fmin, fmax)
        @test length(result.rms) == size(dcs.dyn_ps, 2)
        @test length(result.rms_err) == size(dcs.dyn_ps, 2)
        @test all(result.rms .>= 0)
    end

end
