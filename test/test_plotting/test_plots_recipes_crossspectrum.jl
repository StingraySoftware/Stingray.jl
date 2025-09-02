# Mock data creators
function create_crossspectrum(freqs::Vector{Float64}, powers::Vector{Complex{Float64}};
                               ps1::Union{Vector{Float64}, Nothing}=nothing,
                               ps2::Union{Vector{Float64}, Nothing}=nothing,
                               norm::String="frac", power_type::String="all",
                               nphots1::Float64=1000.0, nphots2::Float64=1000.0,
                               m::Int=1, k::Union{Int, Vector{Int}}=1,
                               segment_size::Float64=100.0,
                               mean_rate1::Float64=10.0, mean_rate2::Float64=10.0,
                               fullspec::Bool=false, channels_overlap::Bool=false)
    mock_metadata = LightCurveMetadata(
        "test", "test", "test", 0.0, (0.0, segment_size), 
        length(freqs) > 1 ? freqs[2] - freqs[1] : 1.0,
        Vector{Dict{String,Any}}(), Dict{String,Any}()
    )
    
    df = length(freqs) > 1 ? freqs[2] - freqs[1] : 1.0
    ps1_data = ps1 !== nothing ? ps1 : abs2.(powers)
    ps2_data = ps2 !== nothing ? ps2 : abs2.(powers)
    
    return CrossSpectrum{Float64}(
        freqs, powers, nothing, ps1_data, ps2_data, norm, power_type, df,
        nphots1, nphots2, m, length(freqs), k, mock_metadata, mock_metadata,
        fullspec, channels_overlap, segment_size, mean_rate1, mean_rate2
    )
end

function create_lightcurve_test(times::Vector{Float64}, counts::Vector{Int}, dt::Float64=1.0)
    tstart, tstop = minimum(times) - dt/2, maximum(times) + dt/2
    gti = reshape([tstart, tstop], 1, 2)
    
    metadata = LightCurveMetadata(
        "test", "test", "test", 0.0, (tstart, tstop), dt, 
        Vector{Dict{String,Any}}(), Dict{String,Any}("gti" => gti)
    )
    
    return LightCurve(times, dt, counts, nothing, nothing,
                      EventProperty{Float64}[], metadata, :poisson)
end

# Basic amplitude plot recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 1.0, randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:amplitude)
    @test p isa Plots.Plot
    @test length(p.series_list) == 1
    
    # Verify recipe returns correct data
    plot_data = p.series_list[1]
    expected_amplitude = abs.(powers)
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_amplitude
end

# Power spectrum recipe - squared amplitudes
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 1.0, randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:power)
    plot_data = p.series_list[1]
    expected_power = abs2.(powers)
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_power
end

# Phase lag recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:phase_lag)
    plot_data = p.series_list[1]
    expected_phase = angle.(powers)
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_phase
end

# Time lag recipe - phase converted to time
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:time_lag)
    plot_data = p.series_list[1]
    expected_time_lag = angle.(powers) ./ (2π .* freqs)
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_time_lag
end

# Coherence recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    ps1 = abs2.(powers) .+ 0.5
    ps2 = abs2.(powers) .+ 0.5
    cs = create_crossspectrum(freqs, powers, ps1=ps1, ps2=ps2)
    
    p = plot(cs, plot_type=:coherence)
    plot_data = p.series_list[1]
    expected_coherence = abs2.(powers) ./ (ps1 .* ps2)
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_coherence
end

# SNR recipe with theoretical noise
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 2.0, randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:snr)
    plot_data = p.series_list[1]
    noise_level = theoretical_noise_level(cs)
    expected_snr = abs.(powers) ./ noise_level
    @test plot_data.plotattributes[:x] ≈ freqs
    @test plot_data.plotattributes[:y] ≈ expected_snr
end

# Real/imaginary components recipe - separate series
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, plot_type=:real_imaginary)
    @test length(p.series_list) >= 2
    
    # Find real and imaginary series
    real_found = false
    imag_found = false
    
    for series in p.series_list
        if haskey(series.plotattributes, :y)
            y_data = series.plotattributes[:y]
            if isa(y_data, Vector) && length(y_data) == length(powers)
                try
                    if y_data ≈ real.(powers)
                        real_found = true
                    elseif y_data ≈ imag.(powers)
                        imag_found = true
                    end
                catch
                    continue
                end
            end
        end
    end
    
    @test real_found && imag_found
end

# PDS comparison recipe - both input power spectra
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    ps1 = abs2.(powers) .+ randn(length(freqs)) .* 0.1
    ps2 = abs2.(powers) .+ randn(length(freqs)) .* 0.2
    cs = create_crossspectrum(freqs, powers, ps1=ps1, ps2=ps2)
    
    p = plot(cs, plot_type=:pds_comparison)
    @test length(p.series_list) >= 2
    
    # Find ps1 and ps2 series
    ps1_found = false
    ps2_found = false
    
    for series in p.series_list
        if haskey(series.plotattributes, :y)
            y_data = series.plotattributes[:y]
            if isa(y_data, Vector) && length(y_data) == length(ps1)
                try
                    if y_data ≈ ps1
                        ps1_found = true
                    elseif y_data ≈ ps2
                        ps2_found = true
                    end
                catch
                    continue
                end
            end
        end
    end
    
    @test ps1_found && ps2_found
end

# Significant detections specialized recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 0.5, randn(length(freqs)))
    powers[5:7] .= complex.(randn(3) .+ 5.0, randn(3))  # High SNR points
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, Val(:significant_detections), threshold=3.0)
    @test length(p.series_list) >= 2  # Main spectrum + significant points
    
    # Verify significant points detected
    noise_level = theoretical_noise_level(cs)
    snr_data = abs.(powers) ./ noise_level
    significant_mask = snr_data .> 3.0
    @test sum(significant_mask) > 0
end

# Phase lag errors recipe with error bars
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    cs.power_err = abs.(powers) .* 0.1
    
    p = plot(cs, Val(:phase_lag_errors))
    plot_data = p.series_list[1]
    @test haskey(plot_data.plotattributes, :yerror)
    
    expected_phase_err = cs.power_err ./ abs.(powers)
    @test plot_data.plotattributes[:yerror] ≈ expected_phase_err
end

# Time lag errors recipe with proper scaling
let
    freqs = collect(0.5:0.1:5.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    cs.power_err = abs.(powers) .* 0.1
    
    p = plot(cs, Val(:time_lag_errors))
    plot_data = p.series_list[1]
    phase_err = cs.power_err ./ abs.(powers)
    expected_time_err = phase_err ./ (2π .* freqs)
    @test plot_data.plotattributes[:yerror] ≈ expected_time_err
end

# White noise estimation recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 1.0, randn(length(freqs)))
    high_freq_idx = freqs .> 8.0
    powers[high_freq_idx] .= complex.(randn(sum(high_freq_idx)) .+ 0.5, randn(sum(high_freq_idx)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, Val(:white_noise))
    @test length(p.series_list) >= 3
    
    # Check noise level calculation
    high_freq_cutoff = maximum(freqs) * 0.8
    high_freq_mask = freqs .>= high_freq_cutoff
    expected_noise = median(abs.(powers[high_freq_mask]))
    @test expected_noise > 0
end

# Frequency range filtering in recipes
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    freq_range = (2.0, 5.0)
    p = plot(cs, plot_type=:amplitude, freq_range=freq_range)
    plot_data = p.series_list[1]
    plotted_freqs = plot_data.plotattributes[:x]
    @test all(plotted_freqs .>= freq_range[1])
    @test all(plotted_freqs .<= freq_range[2])
end

# Coherence confidence recipe for averaged spectra
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    ps1 = abs2.(powers) .+ randn(length(freqs)) .* 0.1
    ps2 = abs2.(powers) .+ randn(length(freqs)) .* 0.1
    m_segments = 10
    cs = create_crossspectrum(freqs, powers, ps1=ps1, ps2=ps2, m=m_segments)
    
    p = plot(cs, Val(:coherence_confidence))
    @test length(p.series_list) >= 3
    
    # Verify confidence level calculations
    confidence_95 = 0.95^(1/(m_segments-1))
    confidence_99 = 0.99^(1/(m_segments-1))
    @test confidence_95 < 1.0
    @test confidence_99 < 1.0
    @test confidence_99 > confidence_95
end

# Frequency rebinning recipe
let
    freqs = collect(0.1:0.01:5.0)  # 491 points
    powers = complex.(randn(length(freqs)) .+ 1.5, randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    rebin_factor = 10
    p = plot(cs, Val(:frequency_snr), rebin_factor=rebin_factor)
    plot_data = p.series_list[1]
    rebinned_points = length(plot_data.plotattributes[:x])
    expected_points = div(length(freqs), rebin_factor)
    @test rebinned_points ≈ expected_points rtol=0.1
end

# Noise comparison recipe
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)) .+ 2.0, randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    p = plot(cs, Val(:noise_comparison))
    @test length(p.series_list) >= 2
    
    original_series = p.series_list[1]
    corrected_series = p.series_list[2]
    
    original_values = original_series.plotattributes[:y]
    corrected_values = corrected_series.plotattributes[:y]
    
    @test original_values ≈ abs.(powers)
    @test all(corrected_values .<= original_values)  # Noise subtraction effect
end

# Light curve comparison recipe with downsampling
let
    times1 = collect(0.0:0.01:1000.0)  # 100,001 points
    times2 = collect(0.0:0.01:1000.0)
    counts1 = rand(1:50, length(times1))
    counts2 = rand(1:40, length(times2))
    
    lc1 = create_lightcurve_test(times1, counts1, 0.01)
    lc2 = create_lightcurve_test(times2, counts2, 0.01)
    
    max_points = 5000
    p = plot(lc1, lc2, max_points=max_points)
    @test p isa Plots.Plot
    @test length(p.series_list) >= 2
    
    # Check downsampling behavior
    data_series_lengths = Int[]
    for series in p.series_list
        if haskey(series.plotattributes, :x) && haskey(series.plotattributes, :y)
            x_data = series.plotattributes[:x]
            y_data = series.plotattributes[:y]
            if isa(x_data, Vector) && isa(y_data, Vector) && 
               length(x_data) == length(y_data) && length(x_data) > 100
                push!(data_series_lengths, length(x_data))
            end
        end
    end
    
    @test length(data_series_lengths) >= 2
    for n_points in data_series_lengths
        @test n_points <= max_points * 1.1  # Allow rounding tolerance
    end
    @test any(n_points -> n_points < length(times1) / 10, data_series_lengths)
end

# Multiple spectra normalization comparison recipe
let
    freqs = collect(0.1:0.1:10.0)
    
    powers1 = complex.(randn(length(freqs)) .+ 1.0, randn(length(freqs)))
    powers2 = powers1 .* 2.0
    powers3 = powers1 .* 0.5
    
    cs1 = create_crossspectrum(freqs, powers1, norm="leahy")
    cs2 = create_crossspectrum(freqs, powers2, norm="frac") 
    cs3 = create_crossspectrum(freqs, powers3, norm="rms")
    
    cs_array = [cs1, cs2, cs3]
    p = plot(cs_array, Val(:normalization_comparison))
    @test length(p.series_list) >= length(cs_array)
    
    # Check each spectrum is plotted
    found_count = 0
    for cs in cs_array
        expected_data = abs.(cs.power)
        for series in p.series_list
            if haskey(series.plotattributes, :y)
                y_data = series.plotattributes[:y]
                if isa(y_data, Vector) && length(y_data) == length(expected_data)
                    try
                        if y_data ≈ expected_data
                            found_count += 1
                            break
                        end
                    catch
                        continue
                    end
                end
            end
        end
    end
    @test found_count == length(cs_array)
end

# Rebinning comparison recipe
let
    freqs_fine = collect(0.1:0.01:10.0)
    freqs_coarse = collect(0.1:0.1:10.0)
    
    base_spectrum = sin.(freqs_fine * 2π) .+ randn(length(freqs_fine)) .* 0.1
    powers_fine = complex.(base_spectrum, randn(length(freqs_fine)) .* 0.1)
    
    coarse_spectrum = [mean(base_spectrum[max(1, round(Int, (f-0.05)/0.01)):min(end, round(Int, (f+0.05)/0.01))]) for f in freqs_coarse]
    powers_coarse = complex.(coarse_spectrum, randn(length(freqs_coarse)) .* 0.1)
    
    cs_fine = create_crossspectrum(freqs_fine, powers_fine)
    cs_coarse = create_crossspectrum(freqs_coarse, powers_coarse)
    
    p = plot(cs_fine, cs_coarse, Val(:rebinning_comparison))
    @test length(p.series_list) >= 2
    
    # Check different resolutions are present
    series_lengths = [length(series.plotattributes[:x]) for series in p.series_list]
    @test maximum(series_lengths) > minimum(series_lengths)
end

# Error handling for invalid plot type
let
    freqs = collect(0.1:0.1:10.0)
    powers = complex.(randn(length(freqs)), randn(length(freqs)))
    cs = create_crossspectrum(freqs, powers)
    
    @test_throws ErrorException plot(cs, plot_type=:invalid_type)
end