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

    @testset "rebin - linear" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata("[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(), test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
        rng = MersenneTwister(42)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = EventList(times1, nothing, meta)
        ev2 = EventList(times2, nothing, meta)
        cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm="frac", silent=true)

        # Test df_new
        cs_rb = rebin(cs, 5 * cs.df; method=:mean)
        @test length(cs_rb.freq) < length(cs.freq)
        @test cs_rb.df ≈ 5 * cs.df
        @test eltype(cs_rb.power) <: Complex
        @test cs_rb.m isa Int
        @test cs_rb.k == 5

        # Test f
        cs_rb2 = rebin(cs, 1.0; f=5.0)
        @test cs_rb2.df ≈ 5 * cs.df
    end

    @testset "rebin_log" begin
        test_gti = [0.0 10.0]
        meta = FITSMetadata("[test]", 1, nothing, Dict{String,Vector}(), Dict{String,Any}(), test_gti, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
        rng = MersenneTwister(42)
        times1 = sort(rand(rng, Uniform(0.0, 10.0), 1000))
        times2 = sort(rand(rng, Uniform(0.0, 10.0), 800))
        ev1 = EventList(times1, nothing, meta)
        ev2 = EventList(times2, nothing, meta)
        cs = Crossspectrum(ev1, ev2, 10.0, 0.001; norm="frac", silent=true)

        cs_log = rebin_log(cs, 0.01)
        @test length(cs_log.freq) < length(cs.freq)
        @test cs_log.k isa Vector{Int}
        @test cs_log.m isa Vector{Int}
        @test eltype(cs_log.power) <: Complex
    end

end
