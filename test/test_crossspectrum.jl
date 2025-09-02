# Helper function to create mock CrossSpectrum for testing
function create_test_crossspectrum(
    freqs::Vector{Float64},
    powers::Vector{Complex{Float64}};
    ps1::Union{Vector{Float64},Nothing} = nothing,
    ps2::Union{Vector{Float64},Nothing} = nothing,
    norm::String = "frac",
    power_type::String = "all",
    nphots1::Float64 = 1000.0,
    nphots2::Float64 = 1000.0,
    m::Int = 1,
    k::Union{Int,Vector{Int}} = 1,
    segment_size::Float64 = 100.0,
    mean_rate1::Float64 = 10.0,
    mean_rate2::Float64 = 10.0,
    fullspec::Bool = false,
    channels_overlap::Bool = false,
)

    mock_metadata = LightCurveMetadata(
        "test",
        "test",
        "test",
        0.0,
        (0.0, segment_size),
        length(freqs) > 1 ? freqs[2] - freqs[1] : 1.0,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}(),
    )

    df = if length(freqs) > 1
        freqs[2] - freqs[1]
    else
        1.0
    end

    ps1_data = ps1 !== nothing ? ps1 : abs2.(powers)
    ps2_data = ps2 !== nothing ? ps2 : abs2.(powers)

    return CrossSpectrum{Float64}(
        freqs,
        powers,
        nothing,
        ps1_data,
        ps2_data,
        norm,
        power_type,
        df,
        nphots1,
        nphots2,
        m,
        length(freqs),
        k,
        mock_metadata,
        mock_metadata,
        fullspec,
        channels_overlap,
        segment_size,
        mean_rate1,
        mean_rate2,
    )
end
# Helper function to create EventList for testing
function create_test_cseventlist(
    times::Vector{Float64},
    energies::Union{Vector{Float64},Nothing} = nothing,
)
    mock_headers =
        FITSIO.FITSHeader(["MJDREF"], [50000.0], ["Modified Julian Date reference"])
    gti = reshape([minimum(times), maximum(times)], 1, 2)

    mock_metadata = FITSMetadata(
        "test.fits",
        1,
        "EVENTS",
        Dict{String,Vector}(),
        mock_headers,
        gti,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

    return EventList(times, energies, mock_metadata)
end

# Helper function to create LightCurve for testing
function create_test_cslightcurve(
    times::Vector{Float64},
    counts::Vector{Int},
    dt::Float64 = 1.0,
)
    tstart, tstop = minimum(times) - dt / 2, maximum(times) + dt / 2
    gti = reshape([tstart, tstop], 1, 2)

    metadata = LightCurveMetadata(
        "test",
        "test",
        "test",
        0.0,
        (tstart, tstop),
        dt,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => gti),
    )

    return LightCurve(
        times,
        dt,
        counts,
        nothing,
        nothing,
        EventProperty{Float64}[],
        metadata,
        :poisson,
    )
end
# Test CrossSpectrum struct construction and helper functions
let
    cs_type = CrossSpectrum{Float64}
    @test cs_type <: AbstractCrossSpectrum{Float64}

    mock_cs_single =
        create_test_crossspectrum([1.0, 2.0, 3.0], [1.0 + 0im, 2.0 + 1im, 3.0 - 1im], m = 1)

    mock_cs_averaged =
        create_test_crossspectrum([1.0, 2.0, 3.0], [1.0 + 0im, 2.0 + 1im, 3.0 - 1im], m = 5)

    @test is_single(mock_cs_single)
    @test !is_averaged(mock_cs_single)
    @test !is_single(mock_cs_averaged)
    @test is_averaged(mock_cs_averaged)
end

# Test theoretical noise level calculation
let
    cs_single = create_test_crossspectrum(
        [1.0, 2.0],
        [1.0 + 0im, 2.0 + 1im],
        nphots1 = 100.0,
        nphots2 = 200.0,
        m = 1,
    )

    cs_averaged = create_test_crossspectrum(
        [1.0, 2.0],
        [1.0 + 0im, 2.0 + 1im],
        nphots1 = 500.0,
        nphots2 = 1000.0,
        m = 5,
    )

    noise_single = theoretical_noise_level(cs_single)
    noise_averaged = theoretical_noise_level(cs_averaged)

    @test noise_single > 0
    @test noise_averaged > 0
    @test noise_averaged < noise_single
end

# Test fill_errors! function
let
    cs =
        create_test_crossspectrum([1.0, 2.0, 3.0], [1.0 + 0im, 2.0 + 1im, 3.0 - 1im], m = 5)

    @test isnothing(cs.power_err)
    fill_errors!(cs)
    @test !isnothing(cs.power_err)
    @test length(cs.power_err) == length(cs.power)
    @test all(cs.power_err .>= 0)
end

# Test white noise level estimation
let
    freqs = collect(0.1:0.1:10.0)
    powers = [complex(1.0, 0.5) for _ in freqs]
    powers[1:20] .*= 3.0

    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, r"Estimated noise level .* much higher than theoretical") begin
        noise_level = white_noise_level(cs)
        @test noise_level > 0
        @test noise_level < mean(abs.(powers[1:20]))
    end
end

# Test noise corrected power
let
    freqs = collect(0.1:0.1:5.0)
    powers = [complex(2.0 + 0.5 * randn(), 0.5 * randn()) for _ in freqs]

    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, r"Estimated noise level .* much higher than theoretical") begin
        corrected = noise_corrected_power(cs)
        @test length(corrected) == length(powers)
        @test all(corrected .>= 0)
    end
end

# Test signal-to-noise ratio calculation
let
    freqs = [1.0, 2.0, 3.0]
    powers = [complex(5.0, 1.0), complex(2.0, 0.5), complex(10.0, 2.0)]

    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, "Too few high frequency points for reliable noise estimation") begin
        snr = signal_to_noise_ratio(cs)
        @test length(snr) == 3
        @test all(snr .>= 0)
        @test snr[3] > snr[2]
    end
end

# Test aliasing detection
let
    freqs = collect(0.1:0.1:5.0)
    powers = [complex(1.0, 0.1) for _ in freqs]

    cs_normal = create_test_crossspectrum(freqs, powers)

    # Test without expecting specific warning logs since detect_aliasing might not produce them
    aliasing_detected, message = detect_aliasing(cs_normal)
    @test !aliasing_detected

    powers_alias = copy(powers)
    high_freq_indices = Int(round(0.8 * length(freqs))):length(freqs)
    powers_alias[high_freq_indices] .*= 20.0

    cs_alias = create_test_crossspectrum(freqs, powers_alias)

    # Test without expecting specific warning logs
    aliasing_detected_2, message_2 = detect_aliasing(cs_alias)
    @test aliasing_detected_2
end

# Test coherence calculation
let
    freqs = [1.0, 2.0, 3.0]
    cross_powers = [complex(2.0, 1.0), complex(1.5, 0.5), complex(3.0, 1.5)]
    ps1 = [4.0, 2.25, 9.0]
    ps2 = [1.0, 1.0, 1.0]

    cs = create_test_crossspectrum(freqs, cross_powers, ps1 = ps1, ps2 = ps2)

    coh = coherence(cs)
    @test length(coh) == 3
    @test all(coh .>= 0.0)
    @test all(coh .<= 1.0)

    cs_avg = create_test_crossspectrum(freqs, cross_powers, ps1 = ps1, ps2 = ps2, m = 5)

    coh_avg = coherence(cs_avg)
    @test length(coh_avg) == 3
    @test all(coh_avg .>= 0.0)
end

# Test phase and time lag calculations
let
    freqs = [1.0, 2.0, 4.0]
    cross_powers = [complex(1.0, 1.0), complex(0.0, 2.0), complex(-1.0, 1.0)]

    cs = create_test_crossspectrum(freqs, cross_powers)

    phases = phase_lag(cs)
    @test length(phases) == 3
    @test all(-π .<= phases .<= π)
    @test phases[1] ≈ π / 4 atol = 1e-10
    @test phases[2] ≈ π / 2 atol = 1e-10

    time_lags = time_lag(cs)
    @test length(time_lags) == 3
    @test time_lags[1] ≈ phases[1] / (2π * freqs[1]) atol = 1e-10
    @test time_lags[2] ≈ phases[2] / (2π * freqs[2]) atol = 1e-10
end

# Test noise properties analysis
let
    freqs = collect(0.5:0.5:10.0)
    powers = [complex(1.0 + 0.1 * randn(), 0.1 * randn()) for _ in freqs]

    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, r"Estimated noise level .* much higher than theoretical") match_mode =
        :any begin
        props = noise_properties(cs)

        required_keys = [
            "theoretical_noise",
            "noise_level_1",
            "noise_level_2",
            "mean_power",
            "std_power",
            "mean_snr",
            "significant_detections",
            "total_frequencies",
            "high_freq_power",
            "noise_to_signal_ratio",
            "mean_rate_1",
            "mean_rate_2",
            "is_averaged",
            "segments_averaged",
            "averaging_improvement",
        ]

        for key in required_keys
            @test haskey(props, key)
        end

        @test props["total_frequencies"] == length(freqs)
        @test props["is_averaged"] == false
        @test props["segments_averaged"] == 1
        @test props["averaging_improvement"] == 1.0
    end
end

# Test significant frequencies detection
let
    freqs = [0.5, 1.0, 1.5, 2.0, 2.5]
    powers = [
        complex(0.1, 0.05),
        complex(5.0, 2.0),
        complex(0.2, 0.1),
        complex(8.0, 3.0),
        complex(0.15, 0.08),
    ]

    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, "Too few high frequency points for reliable noise estimation") match_mode =
        :any begin
        sig_freqs_3 = significant_frequencies(cs, 3.0)
        sig_freqs_5 = significant_frequencies(cs, 5.0)

        @test length(sig_freqs_5) <= length(sig_freqs_3)
    end
end

# Test get_noise_level with different methods
let
    freqs = collect(0.1:0.1:5.0)
    powers = [complex(1.0, 0.1) for _ in freqs]

    cs = create_test_crossspectrum(freqs, powers)

    theoretical_noise = get_noise_level(cs, :theoretical)
    @test theoretical_noise > 0

    @test_logs (:warn, r"Estimated noise level .* much higher than theoretical") begin
        empirical_noise = get_noise_level(cs, :empirical)
        @test empirical_noise > 0
    end

    @test_throws ArgumentError get_noise_level(cs, :invalid_method)
end

# Test quality metrics calculation
let
    freqs = [1.0, 2.0, 3.0, 4.0]
    powers = [complex(5.0, 1.0), complex(2.0, 0.5), complex(8.0, 2.0), complex(1.0, 0.2)]

    cs_single = create_test_crossspectrum(freqs, powers, m = 1)
    cs_averaged = create_test_crossspectrum(freqs, powers, m = 10)

    @test_logs (:warn, "Too few high frequency points for reliable noise estimation") match_mode =
        :any begin
        metrics_single = quality_metrics(cs_single)
        metrics_averaged = quality_metrics(cs_averaged)

        required_keys = [
            "mean_snr",
            "max_snr",
            "significant_fraction",
            "noise_level",
            "dynamic_range",
            "is_averaged",
        ]

        for key in required_keys
            @test haskey(metrics_single, key)
            @test haskey(metrics_averaged, key)
        end

        @test metrics_single["is_averaged"] == false
        @test metrics_averaged["is_averaged"] == true
        @test haskey(metrics_averaged, "segments_averaged")
        @test haskey(metrics_averaged, "expected_improvement")
        @test metrics_averaged["segments_averaged"] == 10
        @test metrics_averaged["expected_improvement"] ≈ sqrt(10) atol = 1e-10
    end
end

# Test rebin function with frequency resolution
let
    freqs = collect(0.5:0.5:10.0)
    powers = [complex(1.0 + 0.1 * i, 0.1 * i) for i = 1:length(freqs)]

    cs = create_test_crossspectrum(freqs, powers)

    rebinned = rebin(cs, 1.0)
    @test rebinned.df ≈ 1.0 atol = 1e-10
    @test length(rebinned.freq) < length(cs.freq)

    @test_throws ErrorException rebin(cs, 0.25)

    rebinned_int = rebin(cs, 2)
    @test rebinned_int.df ≈ 1.0 atol = 1e-10
    @test rebinned.k == 2
end

# Test logarithmic rebinning
let
    freqs = collect(0.1:0.1:10.0)
    powers = [complex(1.0, 0.1) for _ in freqs]

    cs = create_test_crossspectrum(freqs, powers, m = 5)

    log_rebinned = rebin_log(cs; f = 0.1)
    @test length(log_rebinned.freq) < length(cs.freq)
    @test all(log_rebinned.freq .> 0)
    @test issorted(log_rebinned.freq)

    @test_throws ErrorException rebin_log(cs; f = 0.0)
    @test_throws ErrorException rebin_log(cs; f = 1.0)
end

# Test geometric rebinning
let
    freqs = collect(0.1:0.1:5.0)
    powers = [complex(1.0, 0.1) for _ in freqs]

    cs = create_test_crossspectrum(freqs, powers)

    geom_rebinned = geometric_rebin(cs, 1.5)
    @test length(geom_rebinned.freq) < length(cs.freq)
    @test all(geom_rebinned.freq .> 0)
    @test issorted(geom_rebinned.freq)

    @test_throws ErrorException geometric_rebin(cs, 0.8)
end

# Test rebinning helper functions
let
    freqs = [1.0, 2.0, 3.0, 4.0]
    powers = [complex(1.0, 0.1), complex(2.0, 0.2), complex(3.0, 0.3), complex(4.0, 0.4)]

    cs_single = create_test_crossspectrum(freqs, powers, m = 1, k = 3)
    cs_averaged = create_test_crossspectrum(freqs, powers, m = 5, k = [1, 2, 3, 4])

    @test is_rebinned(cs_single)
    @test is_rebinned(cs_averaged)

    samples_single = effective_samples_per_bin(cs_single)
    samples_averaged = effective_samples_per_bin(cs_averaged)

    @test samples_single == 3
    @test samples_averaged == [5, 10, 15, 20]
end

# Test adaptive rebinning
let
    freqs = collect(0.5:0.5:10.0)
    powers = [complex(1.0, 0.1) for _ in freqs]
    cs = create_test_crossspectrum(freqs, powers)

    @test_logs (:warn, r"Estimated noise level .* much higher than theoretical") match_mode =
        :any begin
        adaptive_rebinned = adaptive_rebin(cs, 5.0, 8)
        @test length(adaptive_rebinned.freq) <= length(cs.freq)
        @test length(adaptive_rebinned.freq) > 0

        heavily_rebinned = adaptive_rebin(cs, 100.0, 3)
        @test length(heavily_rebinned.freq) <= length(cs.freq)

        max_rebinned = adaptive_rebin(cs, 1000.0, 5)
        @test length(max_rebinned.freq) < length(cs.freq)

        result1 = adaptive_rebin(cs, 1.0, 2)
        result2 = adaptive_rebin(cs, 1.0, 2)
        @test length(result1.freq) == length(result2.freq)
    end

    minimal_freqs = [1.0, 2.0]
    minimal_powers = [complex(1.0, 0.1), complex(1.0, 0.1)]
    cs_minimal = create_test_crossspectrum(minimal_freqs, minimal_powers)

    @test_logs (:warn, r"Too few high frequency points") begin
        minimal_result = adaptive_rebin(cs_minimal, 3.0, 2)
        @test length(minimal_result.freq) <= length(cs_minimal.freq)
    end
end

# Test Base.show method
let
    cs_single = create_test_crossspectrum([1.0, 2.0], [1.0 + 0im, 2.0 + 1im], m = 1)
    cs_averaged = create_test_crossspectrum([1.0, 2.0], [1.0 + 0im, 2.0 + 1im], m = 5)

    io = IOBuffer()
    show(io, MIME"text/plain"(), cs_single)
    single_output = String(take!(io))
    @test occursin("CrossSpectrum", single_output)
    @test occursin("Single", single_output)

    show(io, MIME"text/plain"(), cs_averaged)
    averaged_output = String(take!(io))
    @test occursin("CrossSpectrum", averaged_output)
    @test occursin("Averaged", averaged_output)
end

# Test CrossSpectrum constructor from EventList
let
    times = collect(0.0:0.1:100.0)
    dt = 0.1

    Random.seed!(42)
    signal1 = 5.0 .* sin.(2π * 0.5 * times) .+ 10.0
    signal2 = 5.0 .* sin.(2π * 0.5 * times .+ π / 6) .+ 10.0

    noise1 = 2.0 .* randn(length(times))
    noise2 = 2.0 .* randn(length(times))

    counts1 = max.(0, round.(Int, signal1 .+ noise1))
    counts2 = max.(0, round.(Int, signal2 .+ noise2))

    event_times1 = Float64[]
    event_times2 = Float64[]

    for i = 1:length(times)
        for _ = 1:counts1[i]
            push!(event_times1, times[i] + rand() * dt - dt / 2)
        end
        for _ = 1:counts2[i]
            push!(event_times2, times[i] + rand() * dt - dt / 2)
        end
    end

    sort!(event_times1)
    sort!(event_times2)

    ev1 = create_test_cseventlist(event_times1)
    ev2 = create_test_cseventlist(event_times2)

    cs = CrossSpectrum(ev1, ev2, 50.0, dt)

    @test is_single(cs)
    @test cs.m == 1
    @test cs.df ≈ 1 / 50.0 atol = 1e-10
    @test cs.segment_size == 50.0
    @test cs.nphots1 > 0
    @test cs.nphots2 > 0
end

# Test AveragedCrossSpectrum from LightCurve
let
    times = collect(0.0:0.1:1000.0)
    dt = 0.1

    Random.seed!(123)
    f1, f2 = 0.5, 2.0

    signal1 = 10.0 .* (sin.(2π * f1 * times) .+ sin.(2π * f2 * times))
    signal2 = 10.0 .* (sin.(2π * f1 * times .+ π / 4) .+ sin.(2π * f2 * times .+ π / 4))

    noise_level = 2.0
    noise1 = noise_level .* randn(length(times))
    noise2 = noise_level .* randn(length(times))

    base_level = 20.0
    counts1 = floor.(Int, (signal1 .+ noise1 .+ base_level))
    counts2 = floor.(Int, (signal2 .+ noise2 .+ base_level))

    counts1 = max.(0, counts1)
    counts2 = max.(0, counts2)

    tstart, tstop = minimum(times), maximum(times)
    gti = [tstart tstop]

    metadata1 = LightCurveMetadata(
        "",
        "",
        "",
        0.0,
        (tstart, tstop),
        dt,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => gti),
    )

    metadata2 = LightCurveMetadata(
        "",
        "",
        "",
        0.0,
        (tstart, tstop),
        dt,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => gti),
    )

    lc1 = LightCurve(
        times,
        dt,
        counts1,
        nothing,
        nothing,
        EventProperty{Float64}[],
        metadata1,
        :poisson,
    )
    lc2 = LightCurve(
        times,
        dt,
        counts2,
        nothing,
        nothing,
        EventProperty{Float64}[],
        metadata2,
        :poisson,
    )

    cs = CrossSpectrum(lc1, lc2, 100.0)
    @test !is_averaged(cs)
    @test is_single(cs)
    @test cs.m == 1

    @test cs.freq[1] ≈ cs.df atol = 1e-10
    @test maximum(cs.freq) ≈ 1 / (2 * dt) atol = 0.1

    theoretical = theoretical_noise_level(cs)
    @test theoretical > 0

    acs = AveragedCrossSpectrum(lc1, lc2, 100.0)
    @test is_averaged(acs)
    @test !is_single(acs)
    @test acs.m > 1

    fill_errors!(acs)
    @test !isnothing(acs.power_err)
    @test all(acs.power_err .>= 0)

    coh = coherence(acs)
    @test all(0 .<= coh .<= 1)

    freq_mask1 = findall(x -> abs(x - f1) < 3 * acs.df, acs.freq)
    freq_mask2 = findall(x -> abs(x - f2) < 3 * acs.df, acs.freq)

    @test !isempty(freq_mask1)
    @test !isempty(freq_mask2)
    @test maximum(coh[freq_mask1]) > 0.2
    @test maximum(coh[freq_mask2]) > 0.2

    phases = phase_lag(acs)
    @test all(-π .<= phases .<= π)

    tlags = time_lag(acs)
    @test length(tlags) == length(acs.freq)

    @test acs.df ≈ 1 / 100.0 atol = 1e-10
    @test length(cs.freq) ≈ length(acs.freq)
    @test cs.df ≈ acs.df atol = 1e-10
end
#all test
let
    # Create time array and simulated data
    times = collect(0.0:0.1:1000.0)
    dt = 0.1

    # Create synthetic signal with better control
    f1, f2 = 0.5, 2.0  # Hz

    # Set random seed for reproducibility
    Random.seed!(123)

    # Create base signals with phase relationship
    signal1 = 10.0 .* (sin.(2π * f1 * times) .+ sin.(2π * f2 * times))
    signal2 = 10.0 .* (sin.(2π * f1 * times .+ π / 4) .+ sin.(2π * f2 * times .+ π / 4))

    # Add controlled noise
    noise_level = 2.0
    noise1 = noise_level .* randn(length(times))
    noise2 = noise_level .* randn(length(times))

    # Ensure positive counts with good SNR
    base_level = 20.0
    counts1 = floor.(Int, (signal1 .+ noise1 .+ base_level))
    counts2 = floor.(Int, (signal2 .+ noise2 .+ base_level))

    # Ensure no negative counts
    counts1 = max.(0, counts1)
    counts2 = max.(0, counts2)

    # Create light curves with proper GTIs
    tstart, tstop = minimum(times), maximum(times)
    gti = [tstart tstop]

    metadata1 = LightCurveMetadata(
        "",
        "",
        "",
        0.0,
        (tstart, tstop),
        dt,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => gti),
    )

    metadata2 = LightCurveMetadata(
        "",
        "",
        "",
        0.0,
        (tstart, tstop),
        dt,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => gti),
    )

    lc1 = LightCurve(
        times,
        dt,
        counts1,
        nothing,
        nothing,
        EventProperty{Float64}[],
        metadata1,
        :poisson,
    )
    lc2 = LightCurve(
        times,
        dt,
        counts2,
        nothing,
        nothing,
        EventProperty{Float64}[],
        metadata2,
        :poisson,
    )

    # Test CrossSpectrum
    cs = CrossSpectrum(lc1, lc2, 100.0)
    @test !is_averaged(cs)
    @test is_single(cs)
    @test cs.m == 1

    # Test frequency properties
    @test cs.freq[1] ≈ cs.df atol = 1e-10
    @test maximum(cs.freq) ≈ 1 / (2 * dt) atol = 0.1

    # Test noise levels
    theoretical = theoretical_noise_level(cs)
    empirical = white_noise_level(cs)
    @test theoretical > 0
    @test empirical > 0

    # Test AveragedCrossSpectrum
    acs = AveragedCrossSpectrum(lc1, lc2, 100.0)
    @test is_averaged(acs)
    @test !is_single(acs)
    @test acs.m > 1

    # Test error estimation
    fill_errors!(acs)
    @test !isnothing(acs.power_err)
    @test all(acs.power_err .>= 0)

    # Test coherence
    coh = coherence(acs)
    @test all(0 .<= coh .<= 1)

    # Find frequency bins near input frequencies
    freq_mask1 = findall(x -> abs(x - f1) < 3 * acs.df, acs.freq)
    freq_mask2 = findall(x -> abs(x - f2) < 3 * acs.df, acs.freq)

    # Test for presence of signal (less stringent coherence test)
    @test !isempty(freq_mask1)
    @test !isempty(freq_mask2)
    @test maximum(coh[freq_mask1]) > 0.2
    @test maximum(coh[freq_mask2]) > 0.2

    # Test significant frequencies detection
    sig_freqs = significant_frequencies(acs)
    @test any(abs.(sig_freqs .- f1) .< 3 * acs.df)
    @test any(abs.(sig_freqs .- f2) .< 3 * acs.df)

    # Test phase lags
    phases = phase_lag(acs)
    @test all(-π .<= phases .<= π)

    # Test time lags
    tlags = time_lag(acs)
    @test length(tlags) == length(acs.freq)

    # Basic property tests
    @test acs.df ≈ 1 / 100.0 atol = 1e-10
    @test length(cs.freq) ≈ length(acs.freq)
    @test cs.df ≈ acs.df atol = 1e-10
end
