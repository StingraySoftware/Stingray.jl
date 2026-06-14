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

        ps = Powerspectrum(lc, 10.0; norm="frac", silent=true)

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

end
