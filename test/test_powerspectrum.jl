using Test

# Test Powerspectrum constructor
function test_powerspectrum()
    dt = 0.1
    times = 0:dt:100
    counts = 100 .+ 10 * sin.(2 * π * 0.5 * times)

    # Test normal powerspectrum calculation
    ps = Powerspectrum(counts, dt)
    @test length(ps.freq) > 0
    @test length(ps.power) > 0
    @test ps.norm == "frac"

    # Test Leahy normalization
    ps_leahy = Powerspectrum(counts, dt, norm="leahy")
    @test length(ps_leahy.freq) > 0
    @test length(ps_leahy.power) > 0
    @test ps_leahy.norm == "leahy"
end

# Test AveragedPowerspectrum constructor
function test_averaged_powerspectrum()
    dt = 0.1
    times = 0:dt:100
    counts = 100 .+ 10 * sin.(2 * π * 0.5 * times)

    # Split counts into segments
    segment_size = 10.0
    n_segments = Int(floor(length(counts) * dt / segment_size))
    segments = [counts[Int((i-1)*segment_size/dt+1):Int(i*segment_size/dt)] for i in 1:n_segments]

    # Test AveragedPowerspectrum calculation
    aps = AveragedPowerspectrum(segments, dt)
    @test length(aps.freq) > 0
    @test length(aps.power) > 0
    @test aps.norm == "frac"

    # Test AveragedPowerspectrum with Leahy normalization
    aps_leahy = AveragedPowerspectrum(segments, dt, norm="leahy")
    @test length(aps_leahy.freq) > 0
    @test length(aps_leahy.power) > 0
    @test aps_leahy.norm == "leahy"
end

# Test DynamicalPowerspectrum constructor
function test_dynamical_powerspectrum()
    dt = 0.1
    times = 0:dt:100
    counts = 100 .+ 10 * sin.(2 * π * 0.5 * times)
    segment_size = 10.0

    # Test DynamicalPowerspectrum calculation
    dps = DynamicalPowerspectrum(counts, dt, segment_size)
    @test length(dps.time) > 0
    @test length(dps.freq) > 0
    @test size(dps.dyn_ps, 1) == length(dps.freq)
    @test size(dps.dyn_ps, 2) == length(dps.time)
    @test dps.norm == "frac"

    # Test DynamicalPowerspectrum with Leahy normalization
    dps_leahy = DynamicalPowerspectrum(counts, dt, segment_size, norm="leahy")
    @test length(dps_leahy.time) > 0
    @test length(dps_leahy.freq) > 0
    @test size(dps_leahy.dyn_ps, 1) == length(dps_leahy.freq)
    @test size(dps_leahy.dyn_ps, 2) == length(dps_leahy.time)
    @test dps_leahy.norm == "leahy"
end

# Run tests
test_powerspectrum()
test_averaged_powerspectrum()
test_dynamical_powerspectrum()

println("All tests passed!")
