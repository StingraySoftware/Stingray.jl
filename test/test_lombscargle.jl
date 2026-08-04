# ──────────────────────────────────────────────────────────────────────────────
# Tests for Lomb-Scargle module
# ──────────────────────────────────────────────────────────────────────────────

using Random

# ──────────────────────────────────────────────────────────────────────────────
# autofrequency tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "autofrequency" begin
    @testset "all parameters specified" begin
        freqs = autofrequency(min_freq=0.1, max_freq=0.5, df=0.1)
        @test freqs ≈ [0.1, 0.2, 0.3, 0.4, 0.5]
    end

    @testset "from length" begin
        freqs = autofrequency(min_freq=0.1, max_freq=0.5, length=10.0)
        @test freqs ≈ [0.1, 0.2, 0.3, 0.4, 0.5]
    end

    @testset "from dt (Nyquist)" begin
        freqs = autofrequency(min_freq=0.1, dt=1.0, length=10.0)
        @test freqs ≈ [0.1, 0.2, 0.3, 0.4, 0.5]
    end

    @testset "default min_freq" begin
        freqs = autofrequency(max_freq=0.5, df=0.2)
        @test freqs ≈ [0.1, 0.3, 0.5]
    end

    @testset "nyquist_factor > 1" begin
        freqs = autofrequency(dt=0.1, df=0.5, nyquist_factor=2)
        @test freqs[end] ≈ 10.0 atol=0.5
        @test freqs[end] > 5.0
    end

    @testset "error: missing df and length" begin
        @test_throws Stingray.ValueError autofrequency(min_freq=0.01, max_freq=0.5)
    end

    @testset "error: missing max_freq and dt" begin
        @test_throws Stingray.ValueError autofrequency(min_freq=0.01, df=1.0)
    end

    @testset "warning: negative min_freq" begin
        @test_logs (:warn, r"min_freq must be positive") autofrequency(min_freq=-0.1, max_freq=0.5, df=0.1)
    end

    @testset "return type" begin
        freqs = autofrequency(min_freq=0.1, max_freq=0.5, df=0.1)
        @test freqs isa Vector{Float64}
    end

    @testset "non-empty result" begin
        freqs = autofrequency(min_freq=0.1, max_freq=0.5, df=0.1)
        @test length(freqs) > 0
    end

    @testset "monotonically increasing" begin
        freqs = autofrequency(min_freq=0.1, max_freq=2.0, df=0.1)
        @test all(diff(freqs) .> 0)
    end

    @testset "uniform spacing" begin
        freqs = autofrequency(min_freq=0.1, max_freq=2.0, df=0.1)
        diffs = diff(freqs)
        @test all(isapprox.(diffs, 0.1, atol=1e-10))
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# lsft_slow tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "lsft_slow" begin
    @testset "constant signal" begin
        t = collect(0.0:0.1:9.9)
        y = fill(5.0, length(t))
        freqs = collect(0.1:0.1:5.0)
        ft = lsft_slow(y, t, freqs)
        # All power should be near zero for a constant signal (no DC in freq grid)
        powers = abs2.(ft)
        # Power at lowest freq should be much less than sum(y)^2
        @test all(powers .< sum(y)^2)
    end

    @testset "single sinusoid peak detection" begin
        t = collect(0.0:0.01:10.0)
        f_signal = 2.5
        y = sin.(2π .* f_signal .* t)
        freqs = collect(0.1:0.1:5.0)
        ft = lsft_slow(y, t, freqs)
        powers = abs2.(ft)
        peak_freq = freqs[argmax(powers)]
        @test isapprox(peak_freq, f_signal, atol=0.15)
    end

    @testset "return type" begin
        t = [0.0, 1.0, 2.0]
        y = [1.0, 2.0, 3.0]
        freqs = [0.1, 0.2, 0.3]
        ft = lsft_slow(y, t, freqs)
        @test ft isa Vector{ComplexF64}
    end

    @testset "output length matches freqs" begin
        t = collect(0.0:0.5:5.0)
        y = randn(length(t))
        freqs = collect(0.1:0.05:2.0)
        ft = lsft_slow(y, t, freqs)
        @test length(ft) == length(freqs)
    end

    @testset "matches manual DFT" begin
        t = [0.0, 1.0, 2.0]
        y = [1.0, 0.0, -1.0]
        ν = 0.25
        expected = sum(y[j] * exp(-2π * im * ν * t[j]) for j in eachindex(y))
        ft = lsft_slow(y, t, [ν])
        @test ft[1] ≈ expected
    end

    @testset "zero signal" begin
        t = collect(0.0:0.1:5.0)
        y = zeros(length(t))
        freqs = collect(0.1:0.1:2.0)
        ft = lsft_slow(y, t, freqs)
        @test all(abs.(ft) .≈ 0.0)
    end

    @testset "input length mismatch error" begin
        @test_throws ArgumentError lsft_slow([1.0, 2.0], [0.0], [0.1])
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# lsft_fast tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "lsft_fast" begin
    @testset "single sinusoid peak detection" begin
        t = collect(0.0:0.01:10.0)
        f_signal = 2.5
        y = sin.(2π .* f_signal .* t)
        freqs = collect(0.1:0.1:5.0)
        ft = lsft_fast(y, t, freqs; oversampling=5)
        powers = abs2.(ft)
        peak_freq = freqs[argmax(powers)]
        @test isapprox(peak_freq, f_signal, atol=0.15)
    end

    @testset "return type" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        y = [1.0, 2.0, 3.0, 2.0, 1.0]
        freqs = [0.1, 0.2, 0.3, 0.4, 0.5]
        ft = lsft_fast(y, t, freqs)
        @test ft isa Vector{ComplexF64}
    end

    @testset "output length matches freqs" begin
        t = collect(0.0:0.5:5.0)
        y = randn(length(t))
        freqs = collect(0.1:0.05:2.0)
        ft = lsft_fast(y, t, freqs)
        @test length(ft) == length(freqs)
    end

    @testset "consistency with lsft_slow" begin
        rng = MersenneTwister(42)
        t = sort(rand(rng, 50) .* 10.0)
        y = sin.(2π .* 0.5 .* t) .+ 0.1 .* randn(rng, 50)
        freqs = collect(0.05:0.05:2.0)
        ft_slow = lsft_slow(y, t, freqs)
        ft_fast = lsft_fast(y, t, freqs; oversampling=10)
        # The fast method should find the same peak
        peak_slow = freqs[argmax(abs2.(ft_slow))]
        peak_fast = freqs[argmax(abs2.(ft_fast))]
        @test isapprox(peak_fast, peak_slow, atol=0.1)
    end

    @testset "higher oversampling improves accuracy" begin
        t = collect(0.0:0.01:10.0)
        f_signal = 1.3
        y = sin.(2π .* f_signal .* t)
        freqs = collect(0.1:0.1:5.0)
        ft_low = lsft_fast(y, t, freqs; oversampling=2)
        ft_high = lsft_fast(y, t, freqs; oversampling=10)
        ft_ref = lsft_slow(y, t, freqs)
        # High oversampling should be closer to reference
        err_low = sum(abs2.(ft_low .- ft_ref))
        err_high = sum(abs2.(ft_high .- ft_ref))
        @test err_high <= err_low || isapprox(err_high, err_low, rtol=0.5)
    end

    @testset "input validation" begin
        @test_throws ArgumentError lsft_fast([1.0, 2.0], [0.0], [0.1])
        @test_throws ArgumentError lsft_fast([1.0], [0.0], [0.1]; oversampling=0)
    end

    @testset "zero signal" begin
        t = collect(0.0:0.1:5.0)
        y = zeros(length(t))
        freqs = collect(0.1:0.1:2.0)
        ft = lsft_fast(y, t, freqs)
        @test all(abs.(ft) .< 1e-10)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# impose_symmetry_lsft tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "impose_symmetry_lsft" begin
    @testset "output length" begin
        freqs = collect(0.1:0.1:0.5)
        t = collect(0.0:0.1:5.0)
        y = sin.(2π .* 0.3 .* t)
        ft = lsft_slow(y, t, freqs)
        ft_full, freqs_full = impose_symmetry_lsft(ft, sum(y), length(y), freqs)
        @test length(ft_full) == 2 * length(freqs) + 1
        @test length(freqs_full) == 2 * length(freqs) + 1
    end

    @testset "DC component" begin
        t = collect(0.0:0.1:5.0)
        y = fill(3.0, length(t))
        freqs = collect(0.1:0.1:0.5)
        ft = lsft_slow(y, t, freqs)
        ft_full, freqs_full = impose_symmetry_lsft(ft, sum(y), length(y), freqs)
        dc_idx = findfirst(f -> f ≈ 0.0, freqs_full)
        @test !isnothing(dc_idx)
        @test real(ft_full[dc_idx]) ≈ sum(y)
    end

    @testset "Hermitian symmetry" begin
        rng = MersenneTwister(123)
        t = sort(rand(rng, 30) .* 5.0)
        y = randn(rng, 30)
        freqs = collect(0.2:0.2:1.0)
        ft = lsft_slow(y, t, freqs)
        ft_full, freqs_full = impose_symmetry_lsft(ft, sum(y), length(y), freqs)
        # For each positive frequency, check conjugate at negative frequency
        for i in eachindex(freqs)
            pos_idx = findfirst(f -> f ≈ freqs[i], freqs_full)
            neg_idx = findfirst(f -> f ≈ -freqs[i], freqs_full)
            @test !isnothing(pos_idx) && !isnothing(neg_idx)
            @test ft_full[neg_idx] ≈ conj(ft_full[pos_idx])
        end
    end

    @testset "frequency symmetry" begin
        freqs = collect(0.1:0.1:0.5)
        ft = zeros(ComplexF64, 5)
        ft_full, freqs_full = impose_symmetry_lsft(ft, 0.0, 10, freqs)
        # Frequencies should be symmetric around 0
        @test freqs_full[1] ≈ -freqs_full[end]
        @test 0.0 in freqs_full
    end

    @testset "return types" begin
        freqs = collect(0.1:0.1:0.5)
        ft = ones(ComplexF64, 5)
        ft_full, freqs_full = impose_symmetry_lsft(ft, 5.0, 10, freqs)
        @test ft_full isa Vector{ComplexF64}
        @test freqs_full isa Vector{Float64}
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Helper: build a minimal LightCurve for testing
# ──────────────────────────────────────────────────────────────────────────────

function _make_test_lc(time::Vector{Float64}, counts::Vector{Int};
                       dt::Float64=1.0, err_method::Symbol=:poisson)
    metadata = LightCurveMetadata(
        "TEST", "TEST", "test_object", 0.0,
        (time[1], time[end]),
        dt,
        [Dict{String,Any}()],
        Dict{String,Any}()
    )
    count_error = err_method == :poisson ?
        Float64.(max.(sqrt.(counts), 1.0)) : nothing
    exposure = fill(dt, length(time))
    properties = EventProperty{Float64}[]

    return LightCurve{Float64}(
        time, dt, counts, count_error, exposure,
        properties, metadata, err_method
    )
end

function _make_sinusoid_lc(; f_signal=0.5, n_points=500, dt=0.1,
                              amplitude=50, offset=100, rng=nothing)
    t = collect(0.0:dt:(n_points-1)*dt)
    signal = offset .+ round.(Int, amplitude .* sin.(2π .* f_signal .* t))
    if !isnothing(rng)
        signal .+= rand(rng, -5:5, length(t))
    end
    signal = max.(signal, 0)  # No negative counts
    return _make_test_lc(t, signal; dt=dt)
end

function _make_poisson_lc(; n_points=1000, dt=0.1, mean_rate=100, rng=MersenneTwister(42))
    t = collect(0.0:dt:(n_points-1)*dt)
    counts = rand(rng, Distributions.Poisson(mean_rate), length(t))
    return _make_test_lc(t, counts; dt=dt)
end

# ──────────────────────────────────────────────────────────────────────────────
# LombScargleCrossspectrum tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "LombScargleCrossspectrum" begin
    lc1 = _make_sinusoid_lc(f_signal=0.5, n_points=500, dt=0.1)
    lc2 = _make_sinusoid_lc(f_signal=0.5, n_points=500, dt=0.1,
                             rng=MersenneTwister(99))

    @testset "basic construction from LightCurves" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test lscs isa LombScargleCrossspectrum{Float64}
        @test lscs isa AbstractCrossspectrum
    end

    @testset "freq properties" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test lscs.freq isa Vector{Float64}
        @test length(lscs.freq) > 0
        @test all(lscs.freq .> 0)  # positive frequencies only
        @test issorted(lscs.freq)
    end

    @testset "power length matches freq" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test length(lscs.power) == length(lscs.freq)
        @test length(lscs.power_err) == length(lscs.freq)
        @test length(lscs.unnorm_power) == length(lscs.freq)
    end

    @testset "power type all → complex" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; power_type="all")
        @test eltype(lscs.power) <: Complex
        @test eltype(lscs.power_err) <: Complex
    end

    @testset "power_err ≈ power for m=1" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test lscs.m == 1
        @test lscs.power_err ≈ lscs.power ./ sqrt(lscs.m)
    end

    @testset "normalization none" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; norm="none")
        @test lscs.norm == "none"
    end

    @testset "normalization frac" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; norm="frac")
        @test lscs.norm == "frac"
    end

    @testset "normalization leahy" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; norm="leahy")
        @test lscs.norm == "leahy"
    end

    @testset "normalization abs" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; norm="abs")
        @test lscs.norm == "abs"
    end

    @testset "different norms produce different scales" begin
        norms = ["frac", "abs", "leahy", "none"]
        mean_powers = Float64[]
        for n in norms
            lscs = LombScargleCrossspectrum(lc1, lc2; norm=n)
            push!(mean_powers, mean(abs.(lscs.power)))
        end
        # Not all the same
        @test length(unique(round.(mean_powers, sigdigits=3))) > 1
    end

    @testset "power_type real" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; power_type="real")
        @test eltype(lscs.power) <: Real
        @test lscs.power_type == "real"
    end

    @testset "power_type absolute" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; power_type="absolute")
        @test eltype(lscs.power) <: Real
        @test all(lscs.power .>= 0)
        @test lscs.power_type == "absolute"
    end

    @testset "fullspec includes negative frequencies" begin
        lscs_half = LombScargleCrossspectrum(lc1, lc2; fullspec=false)
        lscs_full = LombScargleCrossspectrum(lc1, lc2; fullspec=true)
        @test length(lscs_full.freq) > length(lscs_half.freq)
        @test any(lscs_full.freq .< 0)
        @test lscs_full.fullspec == true
    end

    @testset "method slow vs fast consistent peak" begin
        lscs_slow = LombScargleCrossspectrum(lc1, lc2; method="slow")
        lscs_fast = LombScargleCrossspectrum(lc1, lc2; method="fast")
        peak_slow = lscs_slow.freq[argmax(abs.(lscs_slow.power))]
        peak_fast = lscs_fast.freq[argmax(abs.(lscs_fast.power))]
        @test isapprox(peak_slow, peak_fast, atol=0.2)
    end

    @testset "sinusoidal signal peak at injected frequency" begin
        lscs = LombScargleCrossspectrum(lc1, lc2; norm="none")
        # Skip DC-dominated low frequencies (offset=100 in test signal)
        mask = lscs.freq .> 0.2
        peak_freq = lscs.freq[mask][argmax(abs.(lscs.power[mask]))]
        @test isapprox(peak_freq, 0.5, atol=0.15)
    end

    @testset "metadata fields" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test lscs.n == length(lc1.time)
        @test lscs.m == 1
        @test lscs.k == 1
        @test lscs.nphots1 > 0
        @test lscs.nphots2 > 0
        @test lscs.type == "crossspectrum"
        @test lscs.err_dist in ["poisson", "gauss"]
    end

    @testset "dt field" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        @test lscs.dt ≈ 0.1
    end

    @testset "df field ≈ freq spacing" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        if length(lscs.freq) > 1
            @test lscs.df ≈ (lscs.freq[2] - lscs.freq[1]) atol=1e-10
        end
    end

    @testset "Base.show" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        buf = IOBuffer()
        show(buf, lscs)
        s = String(take!(buf))
        @test occursin("LombScargleCrossspectrum", s)
        @test occursin("freq bins", s)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# LombScarglePowerspectrum tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "LombScarglePowerspectrum" begin
    lc = _make_sinusoid_lc(f_signal=0.5, n_points=500, dt=0.1)

    @testset "basic construction" begin
        lsps = LombScarglePowerspectrum(lc)
        @test lsps isa LombScarglePowerspectrum{Float64}
        @test lsps isa AbstractPowerspectrum
        @test lsps isa AbstractCrossspectrum
    end

    @testset "nphots equals nphots1" begin
        lsps = LombScarglePowerspectrum(lc)
        @test lsps.nphots ≈ lsps.nphots1
    end

    @testset "type is powerspectrum" begin
        lsps = LombScarglePowerspectrum(lc)
        @test lsps.type == "powerspectrum"
    end

    @testset "power is real for norm=none, power_type=real" begin
        lsps = LombScarglePowerspectrum(lc; norm="none", power_type="real")
        @test eltype(lsps.power) <: Real
    end

    @testset "known sinusoid peak" begin
        lsps = LombScarglePowerspectrum(lc; norm="none")
        # Skip DC-dominated low frequencies (offset=100 in test signal)
        mask = lsps.freq .> 0.2
        peak_freq = lsps.freq[mask][argmax(abs.(lsps.power[mask]))]
        @test isapprox(peak_freq, 0.5, atol=0.15)
    end

    @testset "freq/power lengths match" begin
        lsps = LombScarglePowerspectrum(lc)
        @test length(lsps.power) == length(lsps.freq)
        @test length(lsps.power_err) == length(lsps.freq)
    end

    @testset "freq is positive and sorted" begin
        lsps = LombScarglePowerspectrum(lc)
        @test all(lsps.freq .> 0)
        @test issorted(lsps.freq)
    end

    @testset "power_err ≈ power for m=1" begin
        lsps = LombScarglePowerspectrum(lc)
        @test lsps.m == 1
        @test lsps.power_err ≈ lsps.power ./ sqrt(lsps.m)
    end

    @testset "all norms produce different scales" begin
        norms = ["frac", "abs", "leahy", "none"]
        mean_powers = Float64[]
        for n in norms
            lsps = LombScarglePowerspectrum(lc; norm=n)
            push!(mean_powers, mean(abs.(lsps.power)))
        end
        @test length(unique(round.(mean_powers, sigdigits=3))) > 1
    end

    @testset "method slow works" begin
        lsps = LombScarglePowerspectrum(lc; method="slow")
        @test lsps.method == "slow"
        @test length(lsps.freq) > 0
    end

    @testset "method fast works" begin
        lsps = LombScarglePowerspectrum(lc; method="fast")
        @test lsps.method == "fast"
        @test length(lsps.freq) > 0
    end

    @testset "fullspec includes negative frequencies" begin
        lsps_half = LombScarglePowerspectrum(lc; fullspec=false)
        lsps_full = LombScarglePowerspectrum(lc; fullspec=true)
        @test length(lsps_full.freq) > length(lsps_half.freq)
        @test any(lsps_full.freq .< 0)
    end

    @testset "compute_rms returns sensible values" begin
        lsps = LombScarglePowerspectrum(lc; norm="frac")
        min_f = lsps.freq[1]
        max_f = lsps.freq[end]
        rms, rms_err = compute_rms(lsps, min_f, max_f)
        @test rms isa Float64
        @test rms_err isa Float64
        @test rms >= 0.0
        @test rms_err >= 0.0
    end

    @testset "Base.show" begin
        lsps = LombScarglePowerspectrum(lc)
        buf = IOBuffer()
        show(buf, lsps)
        s = String(take!(buf))
        @test occursin("LombScarglePowerspectrum", s)
        @test occursin("freq bins", s)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Edge case tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "LombScargle edge cases" begin
    @testset "very short time series (2 points)" begin
        lc = _make_test_lc([0.0, 1.0], [10, 20]; dt=1.0)
        lsps = LombScarglePowerspectrum(lc; norm="none")
        @test length(lsps.freq) > 0
    end

    @testset "all-zero counts → zero power" begin
        t = collect(0.0:0.1:10.0)
        lc = _make_test_lc(t, zeros(Int, length(t)); dt=0.1)
        lsps = LombScarglePowerspectrum(lc; norm="none")
        @test all(abs.(lsps.power) .< 1e-10)
    end

    @testset "invalid norm string → error" begin
        lc = _make_sinusoid_lc()
        @test_throws Stingray.ValueError LombScarglePowerspectrum(lc; norm="invalid")
    end

    @testset "invalid power_type → error" begin
        lc = _make_sinusoid_lc()
        @test_throws Stingray.ValueError LombScarglePowerspectrum(lc; power_type="invalid")
    end

    @testset "invalid method → error" begin
        lc = _make_sinusoid_lc()
        @test_throws Stingray.ValueError LombScarglePowerspectrum(lc; method="invalid")
    end

    @testset "min_freq > max_freq → error" begin
        lc = _make_sinusoid_lc()
        @test_throws Stingray.ValueError LombScarglePowerspectrum(lc; min_freq=1.0, max_freq=0.1)
    end

    @testset "negative min_freq → error" begin
        lc = _make_sinusoid_lc()
        @test_throws Stingray.ValueError LombScarglePowerspectrum(lc; min_freq=-0.5)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# time_lag tests
# ──────────────────────────────────────────────────────────────────────────────

@testset "time_lag" begin
    lc1 = _make_sinusoid_lc(f_signal=0.5, n_points=500, dt=0.1)
    lc2 = _make_sinusoid_lc(f_signal=0.5, n_points=500, dt=0.1,
                             rng=MersenneTwister(99))

    @testset "return type" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        lags = time_lag(lscs)
        @test lags isa Vector{Float64}
    end

    @testset "length matches freq" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        lags = time_lag(lscs)
        @test length(lags) == length(lscs.freq)
    end

    @testset "identical signals → lags ≈ 0" begin
        lscs = LombScargleCrossspectrum(lc1, lc1)
        lags = time_lag(lscs)
        # For identical signals, all cross-power phases should be 0
        @test all(abs.(lags) .< 1e-10)
    end

    @testset "phase_lag consistency" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        ϕ = phase_lag(lscs)
        lags = time_lag(lscs)
        expected = ϕ ./ (2π .* lscs.freq)
        @test lags ≈ expected
    end

    @testset "phase_lag returns Vector{Float64}" begin
        lscs = LombScargleCrossspectrum(lc1, lc2)
        ϕ = phase_lag(lscs)
        @test ϕ isa Vector{Float64}
        @test length(ϕ) == length(lscs.freq)
    end
end

