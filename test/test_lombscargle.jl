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
