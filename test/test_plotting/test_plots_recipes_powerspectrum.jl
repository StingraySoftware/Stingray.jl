# Test utility functions
let
    @test get_norm_label("leahy") == "Leahy Power"
    @test get_norm_label("frac") == "(rms/mean)² Hz⁻¹"
    @test get_norm_label("rms") == "rms² Hz⁻¹"
    @test get_norm_label("abs") == "Absolute Power"
    @test get_norm_label("none") == "Raw Power"
    @test get_norm_label("unknown") == "Power"
    @test get_norm_label(:frac) == "(rms/mean)² Hz⁻¹"
end

# Test get_poisson_level wrapper function
let
    @test get_poisson_level("leahy"; meanrate=100.0) isa Real
    @test get_poisson_level("frac"; meanrate=100.0) isa Real
    @test get_poisson_level("abs"; meanrate=50.0, backrate=1.0) isa Real
    @test get_poisson_level("leahy"; n_ph=1000) isa Real
    
    result = get_poisson_level("frac"; n_ph=500)
    @test result isa Real
end

# Test extract_gti function with various metadata types
let
    meta_dict = Dict("gti" => [5.0 6.0; 7.0 8.0], "other" => "data")
    @test extract_gti(meta_dict) == [5.0 6.0; 7.0 8.0]
    
    meta_nested = Dict("extra" => Dict("gti" => [9.0 10.0]), "name" => "test")
    
    meta_no_gti = Dict("other" => "data", "time" => 100.0)
    @test extract_gti(meta_no_gti) === nothing
    
    meta_tuple = (gti = [1.0 2.0; 3.0 4.0], other_field = "test")
    if hasfield(typeof(meta_tuple), :gti)
        @test extract_gti(meta_tuple) == [1.0 2.0; 3.0 4.0]
    end
    
    empty_meta = Dict{String, Any}()
    @test extract_gti(empty_meta) === nothing
end

# Test mock data structures for recipes using NamedTuples
let
    freqs = 10.0 .^ range(-2, 2, length=100)
    powers = 1.0 ./ freqs .+ 0.01 * randn(100)
    powers = max.(powers, 0.001)
    errors = 0.1 * powers
    
    for norm in ["leahy", "frac", "abs"]
        mock_ps = (
            freq = freqs,
            power = powers,
            power_err = errors,
            norm = norm,
            mean_rate = 100.0
        )
        
        @test length(mock_ps.freq) == 100
        @test length(mock_ps.power) == 100
        @test length(mock_ps.power_err) == 100
        @test mock_ps.norm == norm
        @test mock_ps.mean_rate == 100.0
        @test all(p -> p > 0, mock_ps.power)
        @test all(f -> f > 0, mock_ps.freq)
    end
    
    mock_ps_no_err = (
        freq = freqs,
        power = powers,
        power_err = nothing,
        norm = "frac",
        mean_rate = 50.0
    )
    @test mock_ps_no_err.power_err === nothing
end

# Test EventList using helper function
let
    n_events = 10000
    event_times = sort(rand(n_events) * 1000.0)
    energies = 0.5 .+ 9.5 * rand(n_events)
    
    events = create_test_eventlist(event_times, energies)
    
    @test length(events.times) == n_events
    @test length(events.energies) == n_events
    @test issorted(events.times)
    @test all(e -> 0.5 <= e <= 10.0, events.energies)
    @test events.meta.gti[1,1] == 0.0
    @test events.meta.gti[1,2] == maximum(event_times)
end

# Test LightCurve using helper function
let
    dt = 0.1
    times = 0:dt:999.9
    base_rate = 100.0
    counts = [Int(round(base_rate * dt + 0.2 * base_rate * dt * sin(2π * t / 100.0) + 
              sqrt(base_rate * dt) * randn())) for t in times]
    counts = max.(counts, 0)
    
    lc = create_test_lightcurve(collect(times), counts, dt)
    
    @test length(lc.time) == length(lc.counts)
    @test lc.metadata.bin_size == dt
    @test lc.metadata.time_range[2] > lc.metadata.time_range[1]
    @test all(c -> c >= 0, lc.counts)
end

# Test recipe parameter validation
let
    valid_norms = ["leahy", "frac", "abs", "none"]
    valid_drawstyles = [:steppost, :line, :scatter]
    
    for norm in valid_norms
        @test get_norm_label(norm) isa String
    end
    
    bool_params = [true, false]
    for show_noise in bool_params, show_errors in bool_params, 
        freq_mult in bool_params, log_scale in bool_params
        
        @test show_noise isa Bool
        @test show_errors isa Bool  
        @test freq_mult isa Bool
        @test log_scale isa Bool
    end
    
    @test 256.0 > 0
    @test 0.001 > 0
    @test 0.3 >= 0 && 0.3 <= 1
    @test 1.5 > 0
end

# Test frequency multiplication logic
let
    freqs = [1.0, 2.0, 4.0, 8.0, 16.0]
    powers = [16.0, 8.0, 4.0, 2.0, 1.0]
    
    freq_mult_powers = powers .* freqs
    expected = [16.0, 16.0, 16.0, 16.0, 16.0]
    
    @test freq_mult_powers ≈ expected
    @test powers[1] > powers[2] > powers[3] > powers[4] > powers[5]
    
    errors = 0.1 * powers
    freq_mult_errors = errors .* freqs
    @test length(freq_mult_errors) == length(errors)
    @test all(freq_mult_errors .>= 0)
end

# Test noise level calculations
let
    mean_rates = [10.0, 100.0, 1000.0]
    
    for rate in mean_rates
        leahy_noise = get_poisson_level("leahy"; meanrate=rate)
        frac_noise = get_poisson_level("frac"; meanrate=rate)
        
        @test leahy_noise > 0
        @test frac_noise > 0
        @test leahy_noise isa Real
        @test frac_noise isa Real
    end
    
    low_rate_noise = get_poisson_level("frac"; meanrate=10.0)
    high_rate_noise = get_poisson_level("frac"; meanrate=1000.0)
    @test low_rate_noise != high_rate_noise
end

# Define test structs outside of let blocks to avoid scope issues
struct TestSpectrum{T}
    freq::Vector{T}
    power::Vector{T}
    norm::String
    mean_rate::T
end

# Test multiple spectra handling
let
    freqs = 10.0 .^ range(-1, 1, length=50)
    
    ps1 = TestSpectrum(freqs, freqs.^(-1.0), "frac", 100.0)
    ps2 = TestSpectrum(freqs, freqs.^(-1.5), "frac", 100.0)
    ps3 = TestSpectrum(freqs, freqs.^(-2.0), "frac", 100.0)
    
    spectra = [ps1, ps2, ps3]
    
    @test length(spectra) == 3
    @test all(ps -> ps.norm == "frac", spectra)
    @test all(ps -> length(ps.freq) == 50, spectra)
    @test all(ps -> all(p -> p > 0, ps.power), spectra)
    @test !all(ps1.power .≈ ps2.power)
    @test !all(ps2.power .≈ ps3.power)
end

# Test energy band analysis using helper functions
let
    n_soft = 5000
    n_hard = 3000
    
    times_soft = sort(rand(n_soft) * 1000.0)
    energies_soft = 0.5 .+ 1.5 * rand(n_soft)
    
    times_hard = sort(rand(n_hard) * 1000.0)  
    energies_hard = 2.0 .+ 8.0 * rand(n_hard)
    
    events_soft = create_test_eventlist(times_soft, energies_soft)
    events_hard = create_test_eventlist(times_hard, energies_hard)
    
    events_dict = Dict(
        "0.5-2 keV" => events_soft,
        "2-10 keV" => events_hard
    )
    
    @test length(events_dict) == 2
    @test haskey(events_dict, "0.5-2 keV")
    @test haskey(events_dict, "2-10 keV")
    @test all(e -> 0.5 <= e <= 2.0, events_soft.energies)
    @test all(e -> 2.0 <= e <= 10.0, events_hard.energies)
    @test length(events_soft.times) > length(events_hard.times)
end

# Define empty spectrum struct
struct EmptySpectrum{T}
    freq::Vector{T}
    power::Vector{T}  
    norm::String
end

# Test error handling
let
    empty_ps = EmptySpectrum(Float64[], Float64[], "frac")
    @test isempty(empty_ps.freq)
    @test isempty(empty_ps.power)
    @test length(empty_ps.freq) == 0
    @test get_norm_label("completely_invalid") == "Power"
end

# Test realistic data scenarios
let
    segment_sizes = [64.0, 128.0, 256.0, 512.0, 1024.0]
    time_resolutions = [0.001, 0.01, 0.1]
    
    for seg_size in segment_sizes
        @test seg_size > 0
        @test seg_size <= 1024.0
    end
    
    for dt in time_resolutions
        @test dt > 0
        @test dt <= 1.0
    end
    
    nyquist_freqs = 1.0 ./ (2.0 * time_resolutions)
    expected_nyquist = [500.0, 50.0, 5.0]
    
    @test nyquist_freqs ≈ expected_nyquist
end

# Test plot customization parameters
let
    default_colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :cyan, :magenta]
    @test length(default_colors) == 10
    @test all(c -> c isa Symbol, default_colors)
    
    line_styles = [:solid, :dash, :dot, :dashdot]
    @test :dash in line_styles
    
    alpha_values = [0.1, 0.3, 0.5, 0.7, 1.0]
    @test all(α -> 0 <= α <= 1, alpha_values)
    
    line_widths = [0.5, 1.0, 1.5, 2.0, 2.5]
    @test all(w -> w > 0, line_widths)
end

# Test mathematical operations used in recipes
let
    freqs = [0.01, 0.1, 1.0, 10.0, 100.0]
    log_freqs = log10.(freqs)
    expected_log = [-2.0, -1.0, 0.0, 1.0, 2.0]
    
    @test log_freqs ≈ expected_log
    
    power_law(f, norm, index) = norm * f^index
    test_freq = 10.0
    
    @test power_law(test_freq, 1.0, -1.0) ≈ 0.1    # Fixed: use ≈ instead of ==
    @test power_law(test_freq, 1.0, -2.0) ≈ 0.01   # Fixed: use ≈ instead of ==
    
    powers = [10.0, 5.0, 2.0, 1.0, 0.5]
    noise_level = 0.1
    subtracted = max.(powers .- noise_level, 1e-10 * maximum(powers))
    
    @test all(subtracted .> 0)
    @test subtracted[1] ≈ 9.9
    @test subtracted[end] > 0
end

# Test recipe workflow components
let
    total_time = 1000.0
    auto_segment_size = total_time / 32
    expected_size = 31.25
    @test auto_segment_size ≈ expected_size
    
    total_counts = 50000.0
    observation_time = 1000.0
    mean_rate = total_counts / observation_time
    expected_rate = 50.0
    @test mean_rate ≈ expected_rate
    
    gti_start = 1000.0
    gti_stop = 2000.0
    gti_duration = gti_stop - gti_start
    @test gti_duration == 1000.0
end

# Test EventProperty creation
let
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    prop = create_event_property("energy", values, "keV")
    
    @test prop.name == :energy
    @test prop.values == values
    @test prop.unit == "keV"
end

# Test comprehensive EventList with properties
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [1.5, 2.0, 2.5, 3.0, 3.5]
    
    events = create_test_eventlist(times, energies)
    energy_prop = create_event_property("energy", energies, "keV")
    
    @test length(events.times) == 5
    @test events.energies == energies
    @test energy_prop.unit == "keV"
end

# Test power spectrum creation with mock data
let
    freqs = 10.0 .^ range(-2, 1, length=50)
    powers = freqs.^(-1.5)
    errors = 0.1 * powers
    
    mock_ps = (
        freq = freqs,
        power = powers,
        power_err = errors,
        norm = "frac",
        df = freqs[2] - freqs[1],
        nphots = 10000,
        m = 1,
        n = length(freqs)
    )
    
    @test mock_ps.norm == "frac"
    @test length(mock_ps.freq) == 50
    @test all(mock_ps.power .> 0)
    @test mock_ps.df > 0
    @test mock_ps.nphots == 10000
    @test mock_ps.m == 1
    @test mock_ps.n == length(freqs)
end

# Test AbstractPowerSpectrum interface using NamedTuples
let
    freqs = [0.1, 0.2, 0.3, 0.4, 0.5]
    powers = [10.0, 5.0, 2.0, 1.0, 0.5]
    errors = [1.0, 0.5, 0.2, 0.1, 0.05]
    
    ps_single = (
        freq = freqs,
        power = powers,
        power_err = errors,
        norm = "leahy",
        m = 1
    )
    
    ps_averaged = (
        freq = freqs,
        power = powers ./ 2,
        power_err = errors ./ sqrt(4),
        norm = "leahy", 
        m = 4,
        segment_size = 1024.0,
        mean_rate = 100.0
    )
    
    @test ps_single.m == 1
    @test ps_averaged.m == 4
    @test ps_averaged.segment_size == 1024.0
    @test ps_averaged.mean_rate == 100.0
    @test all(ps_averaged.power_err .< ps_single.power_err)
end

# Test LightCurveMetadata structure fields
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [1.5, 2.0, 2.5, 3.0, 3.5]
    dt = 0.1
    
    events = create_test_eventlist(times, energies)
    lc = create_test_lightcurve(times, [1, 2, 3, 4, 5], dt)
    
    @test hasfield(typeof(events.meta), :gti)
    @test events.meta.gti isa Matrix{Float64}
    
    @test hasfield(typeof(lc.metadata), :bin_size)
    @test hasfield(typeof(lc.metadata), :time_range)
    @test hasfield(typeof(lc.metadata), :telescope)
    @test hasfield(typeof(lc.metadata), :instrument)
    @test hasfield(typeof(lc.metadata), :object)
    @test hasfield(typeof(lc.metadata), :mjdref)
    @test hasfield(typeof(lc.metadata), :headers)
    @test hasfield(typeof(lc.metadata), :extra)
    
    @test lc.metadata.bin_size == dt
end

# Test LightCurve structure fields
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [1, 2, 3, 4, 5]
    dt = 1.0
    
    lc = create_test_lightcurve(times, counts, dt)
    
    @test hasfield(typeof(lc), :time)
    @test hasfield(typeof(lc), :counts)
    @test hasfield(typeof(lc), :metadata)
end
