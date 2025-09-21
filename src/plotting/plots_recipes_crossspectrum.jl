"""
# Usage
```julia
# Basic amplitude plot
plot(cs)
plot(cs, plot_type=:amplitude, xscale=:log10)

plot(cs, plot_type=:phase_lag, show_errors=true, freq_range=(0.1, 10))

# Multiple plots
plot(cs, plot_type=:coherence, show_noise_level=true)
```

# Notes
- Returns (freq_plot, data_plot) for single series or Float64[], Float64[] for multi-series
- Uses @series for multiple plot components
- Handles both averaged and single cross spectra
"""
@recipe function f(
    cs::CrossSpectrum{T};
    plot_type = :amplitude,
    xscale = :log10,
    yscale = :log10,
    show_errors = false,
    freq_range = nothing,
    power_range = nothing,
    linewidth = 2,
    fontsize = 8,
    color = :blue,
    show_noise_level = false,
) where {T}

    grid --> true
    minorgrid --> true
    legend --> :topright
    linewidth --> linewidth
    markersize --> 2
    tickfontsize --> fontsize
    guidefontsize --> fontsize
    legendfontsize --> fontsize

    freq_mask = if !isnothing(freq_range)
        (cs.freq .>= freq_range[1]) .& (cs.freq .<= freq_range[2])
    else
        trues(length(cs.freq))
    end

    freq_plot = cs.freq[freq_mask]

    if plot_type == :amplitude
        if is_averaged(cs)
            title --> "Averaged Cross Spectrum Amplitude ($(cs.m) segments)"
        else
            title --> "Single Cross Spectrum Amplitude"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Amplitude"
        xaxis --> xscale
        yaxis --> yscale

        amplitude_plot = abs.(cs.power[freq_mask])

        if !isnothing(power_range)
            ylims --> power_range
        end

        if show_errors && !isnothing(cs.power_err)
            seriestype --> :scatter
            yerror --> cs.power_err[freq_mask]
            color --> color
            markershape --> :circle
            alpha --> 0.7
            label --> "Cross Spectrum Amplitude"
        else
            seriestype --> :line
            color --> color
            alpha --> 0.8
            label --> "Cross Spectrum Amplitude"
        end

        if show_noise_level
            @series begin
                seriestype --> :hline
                y --> [theoretical_noise_level(cs)]
                color --> :red
                linestyle --> :dash
                alpha --> 0.7
                linewidth --> 1
                label --> "Noise Level"
                Float64[], Float64[]
            end
        end

        return freq_plot, amplitude_plot

    elseif plot_type == :power
        if is_averaged(cs)
            title --> "Averaged Cross Power Spectrum ($(cs.m) segments)"
        else
            title --> "Single Cross Power Spectrum"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Power"
        xaxis --> xscale
        yaxis --> yscale

        power_plot = abs2.(cs.power[freq_mask])

        if !isnothing(power_range)
            ylims --> power_range
        end

        if show_errors && !isnothing(cs.power_err)
            seriestype --> :scatter
            yerror --> (cs.power_err[freq_mask] .* 2 .* abs.(cs.power[freq_mask]))
            color --> color
            markershape --> :circle
            label --> "Cross Power"
        else
            seriestype --> :line
            color --> color
            label --> "Cross Power"
        end

        return freq_plot, power_plot

    elseif plot_type == :phase_lag
        if is_averaged(cs)
            title --> "Averaged Phase Lag Spectrum ($(cs.m) segments)"
        else
            title --> "Single Phase Lag Spectrum"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Phase Lag (radians)"
        xaxis --> xscale
        yaxis --> :identity

        phase_plot = angle.(cs.power[freq_mask])

        seriestype --> :line
        color --> color
        label --> "Phase Lag"

        return freq_plot, phase_plot

    elseif plot_type == :time_lag
        if is_averaged(cs)
            title --> "Averaged Time Lag Spectrum ($(cs.m) segments)"
        else
            title --> "Single Time Lag Spectrum"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Time Lag (s)"
        xaxis --> xscale
        yaxis --> :identity

        time_lag_plot = angle.(cs.power[freq_mask]) ./ (2π .* freq_plot)

        seriestype --> :line
        color --> color
        label --> "Time Lag"

        return freq_plot, time_lag_plot

    elseif plot_type == :raw_coherence
        if is_averaged(cs)
            title --> "Averaged Raw Coherence Function ($(cs.m) segments)"
        else
            title --> "Single Raw Coherence Function"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Raw Coherence"
        xaxis --> xscale
        yaxis --> :identity
        ylims --> (0, 1.1)
    
        # Raw coherence calculation (biased)
        raw_coherence_plot =
            abs2.(cs.power[freq_mask]) ./ (cs.ps1[freq_mask] .* cs.ps2[freq_mask])
    
        seriestype --> :line
        color --> color
        label --> "Raw Coherence"
    
        return freq_plot, raw_coherence_plot
    
    elseif plot_type == :coherence
        if is_averaged(cs)
            title --> "Averaged Coherence Function ($(cs.m) segments)"
        else
            title --> "Single Coherence Function"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Coherence"
        xaxis --> xscale
        yaxis --> :identity
        ylims --> (0, 1.1)
    
        # Use the bias-corrected coherence function
        coherence_values = coherence(cs)
        coherence_plot = coherence_values[freq_mask]
    
        seriestype --> :line
        color --> color
        label --> "Coherence"
    
        return freq_plot, coherence_plot

    elseif plot_type == :snr
        if is_averaged(cs)
            title --> "Averaged Signal-to-Noise Ratio ($(cs.m) segments)"
        else
            title --> "Single Signal-to-Noise Ratio"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "S/N Ratio"
        xaxis --> :log10
        yaxis --> :identity

        # Use empirical noise level (consistent with signal_to_noise_ratio function)
        snr_plot = signal_to_noise_ratio(cs)  # This uses white_noise_level internally

        seriestype --> :line
        color --> :purple
        label --> "S/N Ratio (Empirical)"

        return cs.freq, snr_plot
    
    elseif plot_type == :cross_power_raw
        title --> "Cross Power (Raw Complex Values)"
        xlabel --> "Frequency (Hz)"
        ylabel --> "Cross Power"
        xaxis --> xscale
        yaxis --> yscale

        power_raw = cs.power[freq_mask]

        seriestype --> :line
        color --> color
        label --> "Cross Power (Complex Modulus)"

        return freq_plot, abs.(power_raw)

    elseif plot_type == :real_imaginary
        if is_averaged(cs)
            title --> "Averaged Real and Imaginary Parts ($(cs.m) segments)"
        else
            title --> "Single Real and Imaginary Parts"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Cross Spectrum"
        xaxis --> xscale
        yaxis --> :identity

        @series begin
            seriestype --> :line
            color --> :blue
            linewidth --> linewidth
            label --> "Real Part"
            freq_plot, real.(cs.power[freq_mask])
        end

        @series begin
            seriestype --> :line
            color --> :red
            linewidth --> linewidth
            label --> "Imaginary Part"
            freq_plot, imag.(cs.power[freq_mask])
        end

        return Float64[], Float64[]

    elseif plot_type == :pds_comparison
        if is_averaged(cs)
            title --> "Averaged Power Density Spectra ($(cs.m) segments)"
        else
            title --> "Single Power Density Spectra"
        end
        xlabel --> "Frequency (Hz)"
        ylabel --> "Power"
        xaxis --> xscale
        yaxis --> yscale

        @series begin
            seriestype --> :line
            color --> :blue
            linewidth --> linewidth
            label --> "PDS 1"
            freq_plot, cs.ps1[freq_mask]
        end

        @series begin
            seriestype --> :line
            color --> :red
            linewidth --> linewidth
            label --> "PDS 2"
            freq_plot, cs.ps2[freq_mask]
        end

        return Float64[], Float64[]

    else
        error(
            "Unknown plot_type: $plot_type. Available: :amplitude, :power, :phase_lag, :time_lag, :coherence, :snr, :cross_power_raw, :real_imaginary, :pds_comparison",
        )
    end
end

"""
# Usage
```julia
plot(cs, Val(:significant_detections), threshold=5.0)
```

# Notes
- Highlights frequency bins with S/N > threshold
- Uses median of high frequencies for noise estimation
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:significant_detections};
    threshold = 3.0,
    use_empirical_noise = true,  # New parameter for consistency
) where {T}

    if is_averaged(cs)
        title -->
        "Averaged Cross Spectrum - Significant Detections >$(threshold)σ ($(cs.m) segments)"
    else
        title --> "Single Cross Spectrum - Significant Detections >$(threshold)σ"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Cross Spectrum Amplitude"
    xaxis --> :log10
    yaxis --> :log10
    grid --> true
    legend --> :topright

    # Use consistent noise level estimation
    noise_level = if use_empirical_noise
        white_noise_level(cs)  # Same as signal_to_noise_ratio() uses
    else
        theoretical_noise_level(cs)  # Original behavior
    end
    
    snr_data = abs.(cs.power) ./ noise_level
    significant_mask = snr_data .> threshold

    @series begin
        seriestype --> :line
        color --> :gray
        alpha --> 0.5
        linewidth --> 1
        label --> "All frequencies"
        cs.freq, abs.(cs.power)
    end

    if sum(significant_mask) > 0
        @series begin
            seriestype --> :scatter
            color --> :red
            markersize --> 3
            markershape --> :circle
            label --> "Significant ($(sum(significant_mask)) bins)"
            cs.freq[significant_mask], abs.(cs.power[significant_mask])
        end
    end

    @series begin
        seriestype --> :hline
        y --> [noise_level]
        color --> :black
        linestyle --> :dash
        alpha --> 0.7
        label --> use_empirical_noise ? "Empirical Noise Level" : "Theoretical Noise Level"
        Float64[], Float64[]
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot(cs, Val(:phase_lag_errors), freq_range=(0.1, 10))
```

# Notes
- Error propagation: σ_phase = σ_amplitude / amplitude
- Automatically switches between scatter (with errors) and line plots
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:phase_lag_errors};
    freq_range = nothing,
) where {T}

    if is_averaged(cs)
        title --> "Averaged Phase Lag with Error Bars ($(cs.m) segments)"
    else
        title --> "Single Phase Lag with Error Bars"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Phase Lag (radians)"
    xaxis --> :log10
    yaxis --> :identity
    grid --> true

    freq_mask = if !isnothing(freq_range)
        (cs.freq .>= freq_range[1]) .& (cs.freq .<= freq_range[2])
    else
        trues(length(cs.freq))
    end

    freq_plot = cs.freq[freq_mask]
    phase_plot = angle.(cs.power[freq_mask])

    if !isnothing(cs.power_err)
        amplitude = abs.(cs.power[freq_mask])
        phase_err = cs.power_err[freq_mask] ./ amplitude

        seriestype --> :scatter
        yerror --> phase_err
        color --> :blue
        markersize --> 3
        markershape --> :circle
        label --> "Phase Lag with errors"
    else
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Phase Lag"
    end

    if !isnothing(freq_range)
        xlims --> freq_range
    end

    @series begin
        seriestype --> :hline
        y --> [0.0]
        color --> :black
        linestyle --> :dash
        linewidth --> 1
        alpha --> 0.5
        label --> ""
        Float64[], Float64[]
    end

    return freq_plot, phase_plot
end
"""
# Usage
```julia
plot(cs, Val(:time_lag_errors), freq_range=(0.1, 10))
```

# Notes
- Converts phase to time: τ = φ/(2πf)
- Shows error bars when available
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:time_lag_errors};
    freq_range = nothing,
) where {T}

    if is_averaged(cs)
        title --> "Averaged Time Lag with Error Bars ($(cs.m) segments)"
    else
        title --> "Single Time Lag with Error Bars"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Time Lag (s)"
    xaxis --> :log10
    yaxis --> :identity
    grid --> true

    freq_mask = if !isnothing(freq_range)
        (cs.freq .>= freq_range[1]) .& (cs.freq .<= freq_range[2])
    else
        trues(length(cs.freq))
    end

    freq_plot = cs.freq[freq_mask]
    time_lag_plot = angle.(cs.power[freq_mask]) ./ (2π .* freq_plot)

    if !isnothing(cs.power_err)
        amplitude = abs.(cs.power[freq_mask])
        phase_err = cs.power_err[freq_mask] ./ amplitude
        time_err = phase_err ./ (2π .* freq_plot)

        seriestype --> :scatter
        yerror --> time_err
        color --> :blue
        markersize --> 3
        markershape --> :circle
        label --> "Time Lag with errors"
    else
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Time Lag"
    end

    if !isnothing(freq_range)
        xlims --> freq_range
    end

    @series begin
        seriestype --> :hline
        y --> [0.0]
        color --> :black
        linestyle --> :dash
        linewidth --> 1
        alpha --> 0.5
        label --> ""
        Float64[], Float64[]
    end

    return freq_plot, time_lag_plot
end

"""
# Usage
```julia
plot(cs, Val(:white_noise), high_freq_fraction=0.2)
```

# Notes
- Uses median of highest 20% frequencies as noise floor estimate
- Shows noise region and estimated level
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:white_noise};
    high_freq_fraction = 0.2,
) where {T}

    if is_averaged(cs)
        title --> "Averaged Cross Spectrum - White Noise Estimation ($(cs.m) segments)"
    else
        title --> "Single Cross Spectrum - White Noise Estimation"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Power"
    xaxis --> :log10
    yaxis --> :log10
    grid --> true
    legend --> :topright

    freq_cutoff = maximum(cs.freq) * (1.0 - high_freq_fraction)
    high_freq_mask = cs.freq .>= freq_cutoff
    white_noise = median(abs.(cs.power[high_freq_mask]))

    @series begin
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Cross Spectrum"
        cs.freq, abs.(cs.power)
    end

    @series begin
        seriestype --> :vline
        x --> [freq_cutoff]
        color --> :green
        linestyle --> :dash
        alpha --> 0.7
        label --> "Noise Region"
        Float64[], Float64[]
    end

    @series begin
        seriestype --> :hline
        y --> [white_noise]
        color --> :red
        linestyle --> :dash
        linewidth --> 2
        label --> "White Noise Level"
        Float64[], Float64[]
    end

    return Float64[], Float64[]
end

"""
# Usage
```julia
plot(cs, Val(:noise_diagnosis))
```

# Notes
- Comprehensive 2x2 diagnostic panel with aliasing detection
- Shows cross spectrum, S/N ratio, power distribution, and frequency comparison
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:noise_diagnosis}) where {T}

    layout --> (2, 2)
    size --> (1000, 800)

    nyquist_freq = maximum(cs.freq)
    high_freq_mask = cs.freq .> 0.8 * nyquist_freq
    aliasing_detected =
        length(high_freq_mask) > 0 &&
        median(abs.(cs.power[high_freq_mask])) > 5.0 * theoretical_noise_level(cs)

    @series begin
        subplot --> 1
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        alpha --> 0.8
        label --> "Cross Spectrum"
        title --> "Noise Level Comparison"
        xlabel --> "Frequency (Hz)"
        ylabel --> "Power"
        xaxis --> :log10
        yaxis --> :log10
        cs.freq, abs.(cs.power)
    end

    theoretical_noise = theoretical_noise_level(cs)
    high_freq_idx = cs.freq .> maximum(cs.freq) * 0.8
    empirical_noise =
        length(high_freq_idx) > 0 ? median(abs.(cs.power[high_freq_idx])) :
        theoretical_noise

    @series begin
        subplot --> 1
        seriestype --> :hline
        y --> [theoretical_noise]
        color --> :red
        linestyle --> :dash
        linewidth --> 2
        alpha --> 0.7
        label --> "Theoretical Noise"
        Float64[], Float64[]
    end

    @series begin
        subplot --> 1
        seriestype --> :hline
        y --> [empirical_noise]
        color --> :green
        linestyle --> :dot
        linewidth --> 2
        alpha --> 0.7
        label --> "Empirical Noise"
        Float64[], Float64[]
    end

    snr_data = abs.(cs.power) ./ theoretical_noise
    @series begin
        subplot --> 2
        seriestype --> :line
        color --> :purple
        linewidth --> 2
        label --> "S/N Ratio"
        title --> "Signal-to-Noise vs Frequency"
        xlabel --> "Frequency (Hz)"
        ylabel --> "S/N"
        xaxis --> :log10
        yaxis --> :log10
        cs.freq, max.(snr_data, 0.01)
    end

    @series begin
        subplot --> 3
        seriestype --> :histogram
        bins --> 50
        color --> :orange
        alpha --> 0.7
        label --> ""
        title --> "Power Distribution"
        xlabel --> "Log₁₀(Power)"
        ylabel --> "Frequency"
        log10.(abs.(cs.power))
    end

    high_freq_mask = cs.freq .> 0.8 * nyquist_freq
    low_freq_mask = cs.freq .< 0.2 * nyquist_freq

    plot_color = aliasing_detected ? :red : :green
    plot_title = aliasing_detected ? "⚠️ Aliasing Suspected" : "✓ No Obvious Aliasing"

    @series begin
        subplot --> 4
        seriestype --> :scatter
        color --> plot_color
        alpha --> 0.6
        markersize --> 3
        label --> "High Freq (>0.8 Nyquist)"
        title --> plot_title
        xlabel --> "Frequency (Hz)"
        ylabel --> "Power"
        xaxis --> :log10
        yaxis --> :log10
        cs.freq[high_freq_mask], abs.(cs.power[high_freq_mask])
    end

    @series begin
        subplot --> 4
        seriestype --> :scatter
        color --> :blue
        alpha --> 0.6
        markersize --> 3
        label --> "Low Freq (<0.2 Nyquist)"
        cs.freq[low_freq_mask], abs.(cs.power[low_freq_mask])
    end

    return Float64[], Float64[]
end

"""
# Usage
```julia
plot(cs, Val(:analysis))
```

# Notes
- Comprehensive 2x3 analysis dashboard
- Shows amplitude, phase lag, time lag, coherence, real/imaginary parts, and input spectra
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:analysis}) where {T}
    layout --> (2, 3)
    size --> (1200, 800)

    @series begin
        subplot --> 1
        plot_type --> :amplitude
        if is_averaged(cs)
            title --> "Averaged Amplitude ($(cs.m) segments)"
        else
            title --> "Single Amplitude"
        end
        cs
    end

    @series begin
        subplot --> 2
        plot_type --> :phase_lag
        if is_averaged(cs)
            title --> "Averaged Phase Lag ($(cs.m) segments)"
        else
            title --> "Single Phase Lag"
        end
        cs
    end

    @series begin
        subplot --> 3
        plot_type --> :time_lag
        if is_averaged(cs)
            title --> "Averaged Time Lag ($(cs.m) segments)"
        else
            title --> "Single Time Lag"
        end
        cs
    end

    @series begin
        subplot --> 4
        plot_type --> :coherence
        if is_averaged(cs)
            title --> "Averaged Coherence ($(cs.m) segments)"
        else
            title --> "Single Coherence"
        end
        cs
    end

    @series begin
        subplot --> 5
        plot_type --> :real_imaginary
        if is_averaged(cs)
            title --> "Averaged Real & Imaginary ($(cs.m) segments)"
        else
            title --> "Single Real & Imaginary"
        end
        cs
    end

    @series begin
        subplot --> 6
        plot_type --> :pds_comparison
        if is_averaged(cs)
            title --> "Averaged Input Spectra ($(cs.m) segments)"
        else
            title --> "Single Input Spectra"
        end
        cs
    end

    return Float64[], Float64[]
end

"""
# Usage
```julia
plot(cs, Val(:noise_analysis))
```

# Notes
- Detailed 2x2 noise analysis with noise subtraction effects
- Shows original vs noise-corrected power
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:noise_analysis}) where {T}

    layout --> (2, 2)
    size --> (1000, 800)

    noise_level = theoretical_noise_level(cs)
    snr_data = abs.(cs.power) ./ noise_level
    noise_corrected = max.(abs.(cs.power) .- noise_level, 0.01 * noise_level)

    spectrum_info = is_averaged(cs) ? " ($(cs.m) segments)" : ""

    @series begin
        subplot --> 1
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Cross Spectrum"
        title --> "Cross Spectrum vs Noise Level$spectrum_info"
        xlabel --> "Frequency (Hz)"
        ylabel --> "Power"
        xaxis --> :log10
        yaxis --> :log10
        cs.freq, abs.(cs.power)
    end

    @series begin
        subplot --> 1
        seriestype --> :hline
        y --> [noise_level]
        color --> :red
        linestyle --> :dash
        linewidth --> 2
        label --> "Theoretical Noise"
        Float64[], Float64[]
    end

    @series begin
        subplot --> 2
        seriestype --> :line
        color --> :green
        linewidth --> 2
        label --> "S/N Ratio"
        title --> "Signal-to-Noise Ratio$spectrum_info"
        xlabel --> "Frequency (Hz)"
        ylabel --> "S/N"
        xaxis --> :log10
        cs.freq, snr_data
    end

    @series begin
        subplot --> 2
        seriestype --> :hline
        y --> [3.0]
        color --> :black
        linestyle --> :dash
        alpha --> 0.7
        label --> "3σ threshold"
        Float64[], Float64[]
    end

    @series begin
        subplot --> 3
        seriestype --> :line
        color --> :purple
        linewidth --> 2
        label --> "Noise Subtracted"
        title --> "Noise-Corrected Power$spectrum_info"
        xlabel --> "Frequency (Hz)"
        ylabel --> "Corrected Power"
        xaxis --> :log10
        yaxis --> :log10
        cs.freq, noise_corrected
    end

    mean_snr = mean(snr_data)
    significant_bins = sum(snr_data .> 3.0)

    @series begin
        subplot --> 4
        seriestype --> :bar
        color --> [:blue, :red, :green]
        label --> ""
        title --> "Summary Stats$spectrum_info"
        xlabel --> "Metric"
        ylabel --> "Value"
        [1, 2, 3], [mean_snr, noise_level, significant_bins]
    end

    return Float64[], Float64[]
end

"""
# Usage
```julia
plot(cs, Val(:coherence_confidence))
```

# Notes
- Confidence levels: C_α = α^(1/(m-1)) where m is number of segments
- Only meaningful for averaged cross spectra (m > 1)
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:coherence_confidence}) where {T}

    if is_averaged(cs)
        title --> "Averaged Coherence with Confidence Levels ($(cs.m) segments)"
    else
        title --> "Single Coherence with Confidence Levels"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Coherence"
    xaxis --> :log10
    ylims --> (0, 1.1)
    grid --> true

    coh_values = abs2.(cs.power) ./ (cs.ps1 .* cs.ps2)

    m_eff = is_averaged(cs) ? cs.m : 1
    confidence_95 = fill(0.95^(1 / (m_eff - 1)), length(cs.freq))
    confidence_99 = fill(0.99^(1 / (m_eff - 1)), length(cs.freq))

    @series begin
        seriestype --> :line
        color --> :green
        linewidth --> 2
        label --> "Coherence"
        cs.freq, coh_values
    end

    if !isnothing(cs.power_err)
        coh_err = cs.power_err ./ abs.(cs.power)
        @series begin
            seriestype --> :scatter
            yerror --> coh_err
            color --> :green
            alpha --> 0.7
            markersize --> 2
            label --> ""
            cs.freq, coh_values
        end
    end

    @series begin
        seriestype --> :line
        color --> :gray
        linestyle --> :dash
        alpha --> 0.7
        label --> "95% confidence"
        cs.freq, confidence_95
    end

    @series begin
        seriestype --> :line
        color --> :darkgray
        linestyle --> :dot
        alpha --> 0.7
        label --> "99% confidence"
        cs.freq, confidence_99
    end

    return Float64[], Float64[]
end

"""
# Usage
```julia
plot(cs, Val(:frequency_snr), rebin_factor=10)
```

# Notes
- Frequency-resolved S/N with optional rebinning for noise reduction
- Shows 3σ detection threshold
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:frequency_snr}; rebin_factor = 10) where {T}

    if is_averaged(cs)
        title --> "Averaged Frequency-Resolved Signal-to-Noise ($(cs.m) segments)"
    else
        title --> "Single Frequency-Resolved Signal-to-Noise"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "S/N Ratio"
    xaxis --> :log10
    grid --> true

    n_bins = div(length(cs.freq), rebin_factor)
    freq_rebinned = Vector{Float64}(undef, n_bins)
    snr_rebinned = Vector{Float64}(undef, n_bins)

    noise_level = theoretical_noise_level(cs)
    snr_full = abs.(cs.power) ./ noise_level

    for i = 1:n_bins
        start_idx = (i - 1) * rebin_factor + 1
        end_idx = min(i * rebin_factor, length(cs.freq))

        freq_rebinned[i] = mean(cs.freq[start_idx:end_idx])
        snr_rebinned[i] = mean(snr_full[start_idx:end_idx])
    end

    @series begin
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        marker --> :circle
        markersize --> 3
        label --> "Rebinned S/N (factor $rebin_factor)"
        freq_rebinned, snr_rebinned
    end

    @series begin
        seriestype --> :hline
        y --> [1.0]
        color --> :black
        linestyle --> :solid
        alpha --> 0.5
        label --> "Unity"
        Float64[], Float64[]
    end

    @series begin
        seriestype --> :hline
        y --> [3.0]
        color --> :red
        linestyle --> :dash
        alpha --> 0.7
        label --> "3σ detection"
        Float64[], Float64[]
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot(cs, Val(:noise_comparison))
```

# Notes
- Compares original cross spectrum with noise-subtracted version
- Shows theoretical noise level
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:noise_comparison}) where {T}

    if is_averaged(cs)
        title --> "Averaged Cross Spectrum: Original vs Noise-Corrected ($(cs.m) segments)"
    else
        title --> "Single Cross Spectrum: Original vs Noise-Corrected"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Power"
    xaxis --> :log10
    yaxis --> :log10
    grid --> true
    legend --> :topright

    noise_level = theoretical_noise_level(cs)
    noise_corrected = max.(abs.(cs.power) .- noise_level, 0.01 * noise_level)

    @series begin
        seriestype --> :line
        color --> :blue
        alpha --> 0.7
        linewidth --> 2
        label --> "Original"
        cs.freq, abs.(cs.power)
    end

    @series begin
        seriestype --> :line
        color --> :red
        linewidth --> 2
        label --> "Noise Subtracted"
        cs.freq, noise_corrected
    end

    @series begin
        seriestype --> :hline
        y --> [noise_level]
        color --> :gray
        linestyle --> :dash
        alpha --> 0.8
        label --> "Noise Level"
        Float64[], Float64[]
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot([cs1, cs2, cs3], Val(:normalization_comparison), 
     norm_labels=["Leahy", "Fractional RMS", "Absolute RMS"])
```

# Notes
- Compares different normalization schemes
- Layout: stacked plots for each normalization
"""
@recipe function f(
    cs_array::Vector{CrossSpectrum{T}},
    ::Val{:normalization_comparison};
    norm_labels = ["Leahy", "Fractional RMS", "Absolute RMS"],
) where {T}

    layout --> (length(cs_array), 1)
    size --> (800, 300 * length(cs_array))

    for (i, cs) in enumerate(cs_array)
        @series begin
            subplot --> i
            seriestype --> :line
            color --> :black
            linewidth --> 2
            label --> ""
            title --> i <= length(norm_labels) ? norm_labels[i] : "Spectrum $i"
            xlabel --> "Frequency (Hz)"
            ylabel --> i == 1 ? "Leahy cross-power" :
            (i == 2 ? "RMS cross-power" : "Absolute cross-power")
            xaxis --> :identity
            yaxis --> :log10
            cs.freq, abs.(cs.power)
        end
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot(cs, Val(:lag_frequency_errors), 
     freq_range=(0.1, 10), reference_lag=(frequency=1.0, expected_lag=0.1))
```

# Notes
- Shows time lag vs frequency with error bars
- Optional reference lines for expected values
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:lag_frequency_errors};
    rebin_factor = 1,
    freq_range = nothing,
    reference_lag = nothing,
) where {T}

    if is_averaged(cs)
        title --> "Averaged Lag-Frequency Spectrum ($(cs.m) segments)"
    else
        title --> "Single Lag-Frequency Spectrum"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Time Lag (s)"
    xaxis --> :identity
    grid --> true

    freq_mask = if !isnothing(freq_range)
        (cs.freq .>= freq_range[1]) .& (cs.freq .<= freq_range[2])
    else
        trues(length(cs.freq))
    end

    freq_plot = cs.freq[freq_mask]
    lag_plot = angle.(cs.power[freq_mask]) ./ (2π .* freq_plot)

    if !isnothing(freq_range)
        xlims --> freq_range
    end

    @series begin
        seriestype --> :hline
        y --> [0.0]
        color --> :black
        linestyle --> :dash
        linewidth --> 2
        label --> ""
        Float64[], Float64[]
    end

    if !isnothing(reference_lag) && haskey(reference_lag, :frequency)
        @series begin
            seriestype --> :vline
            x --> [reference_lag[:frequency]]
            color --> :gray
            linestyle --> :dash
            alpha --> 0.7
            label --> ""
            Float64[], Float64[]
        end
    end

    if !isnothing(reference_lag) && haskey(reference_lag, :expected_lag)
        @series begin
            seriestype --> :line
            color --> :green
            linewidth --> 2
            label --> "Expected Time Lag"
            freq_plot, fill(reference_lag[:expected_lag], length(freq_plot))
        end
    end

    if !isnothing(cs.power_err)
        amplitude = abs.(cs.power[freq_mask])
        phase_err = cs.power_err[freq_mask] ./ amplitude
        lag_err = phase_err ./ (2π .* freq_plot)

        seriestype --> :scatter
        yerror --> lag_err
        color --> :blue
        markersize --> 3
        markershape --> :circle
        label --> "Measured Time Lag"

        return freq_plot, lag_plot
    else
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Time Lag"

        return freq_plot, lag_plot
    end
end
"""
# Usage
```julia
plot(cs, Val(:phase_frequency_errors), 
     freq_range=(0.1, 10), reference_phase=(frequency=1.0, expected_phase=0.5))
```

# Notes
- Shows phase lag vs frequency with error bars
- Optional reference lines for expected values
"""
@recipe function f(
    cs::CrossSpectrum{T},
    ::Val{:phase_frequency_errors};
    freq_range = nothing,
    reference_phase = nothing,
) where {T}

    if is_averaged(cs)
        title --> "Averaged Phase-Frequency Spectrum ($(cs.m) segments)"
    else
        title --> "Single Phase-Frequency Spectrum"
    end
    xlabel --> "Frequency (Hz)"
    ylabel --> "Phase Lag (rad)"
    xaxis --> :identity
    grid --> true

    freq_mask = if !isnothing(freq_range)
        (cs.freq .>= freq_range[1]) .& (cs.freq .<= freq_range[2])
    else
        trues(length(cs.freq))
    end

    freq_plot = cs.freq[freq_mask]
    phase_plot = angle.(cs.power[freq_mask])

    if !isnothing(freq_range)
        xlims --> freq_range
    end

    @series begin
        seriestype --> :hline
        y --> [0.0]
        color --> :black
        linestyle --> :dash
        linewidth --> 2
        label --> ""
        Float64[], Float64[]
    end

    if !isnothing(reference_phase) && haskey(reference_phase, :frequency)
        @series begin
            seriestype --> :vline
            x --> [reference_phase[:frequency]]
            color --> :gray
            linestyle --> :dash
            alpha --> 0.7
            label --> ""
            Float64[], Float64[]
        end
    end

    if !isnothing(reference_phase) && haskey(reference_phase, :expected_phase)
        @series begin
            seriestype --> :hline
            y --> [reference_phase[:expected_phase]]
            color --> :green
            linewidth --> 2
            label --> "Expected Phase Lag"
            Float64[], Float64[]
        end
    end

    if !isnothing(cs.power_err)
        amplitude = abs.(cs.power[freq_mask])
        phase_err = cs.power_err[freq_mask] ./ amplitude

        seriestype --> :scatter
        yerror --> phase_err
        color --> :blue
        markersize --> 3
        markershape --> :circle
        label --> "Measured Phase Lag"

        return freq_plot, phase_plot
    else
        seriestype --> :line
        color --> :blue
        linewidth --> 2
        label --> "Phase Lag"

        return freq_plot, phase_plot
    end
end

"""
# Usage
```julia
plot(cs_array, energy_bands, Val(:lag_energy))
```

# Notes
- Common in X-ray timing analysis
- Shows energy-dependent time lags
"""
@recipe function f(
    cs_array::Vector{CrossSpectrum{T}},
    energy_bands::Vector{<:Real},
    ::Val{:lag_energy},
) where {T}

    title --> "Lag-Energy Spectrum"
    xlabel --> "Energy (keV)"
    ylabel --> "Time Lag (s)"
    grid --> true

    avg_lags = Float64[]
    lag_errors = Float64[]

    for cs in cs_array
        time_lags = angle.(cs.power) ./ (2π .* cs.freq)
        avg_lag = mean(time_lags)
        lag_error = std(time_lags) / sqrt(length(cs.freq))

        push!(avg_lags, avg_lag)
        push!(lag_errors, lag_error)
    end

    seriestype --> :scatter
    yerror --> lag_errors
    color --> :purple
    markersize --> 5
    markershape --> :circle
    label --> "Energy-dependent Lag"

    return energy_bands[1:end-1], avg_lags
end
"""
# Usage
```julia
plot(cs_original, cs_rebinned, Val(:rebinning_comparison), rebin_type="Linear")
```

# Notes
- Demonstrates effect of rebinning on cross spectrum
- Shows original vs rebinned data
"""
@recipe function f(
    cs_original::CrossSpectrum{T},
    cs_rebinned::CrossSpectrum{T},
    ::Val{:rebinning_comparison};
    rebin_type = "Linear",
) where {T}

    title --> "$rebin_type Rebinning Comparison"
    xlabel --> "Frequency (Hz)"
    ylabel --> "Cross Spectrum Amplitude"
    xaxis --> :log10
    yaxis --> :log10
    grid --> true
    legend --> :topright

    @series begin
        seriestype --> :line
        color --> :lightblue
        alpha --> 0.6
        linewidth --> 1
        label --> "Original"
        cs_original.freq, abs.(cs_original.power)
    end

    @series begin
        seriestype --> :line
        color --> :darkblue
        linewidth --> 2
        label --> "$rebin_type Rebinned"
        cs_rebinned.freq, abs.(cs_rebinned.power)
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot(lc1, lc2, colors=[:blue, :red], max_points=10000, time_range=(0, 100))
```

# Notes
- Automatically downsamples for performance if n > max_points
- Shows point count in legend
"""
@recipe function f(
    lc1::LightCurve,
    lc2::LightCurve;
    colors = [:blue, :red],
    labels = ["Light Curve 1", "Light Curve 2"],
    max_points = 10000,
    time_range = nothing,
)

    title --> "Light Curves Comparison"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    legend --> :topright
    linewidth --> 2

    # Downsample function
    function quick_downsample(time_vec, count_vec, max_pts)
        n = length(time_vec)
        if n <= max_pts
            return time_vec, count_vec
        end
        step = max(1, div(n, max_pts))
        return time_vec[1:step:end], count_vec[1:step:end]
    end

    # Process data
    t1, c1 = quick_downsample(lc1.time, lc1.counts, max_points)
    t2, c2 = quick_downsample(lc2.time, lc2.counts, max_points)

    # Apply time filter
    if !isnothing(time_range)
        xlims --> time_range
        mask1 = (t1 .>= time_range[1]) .& (t1 .<= time_range[2])
        mask2 = (t2 .>= time_range[1]) .& (t2 .<= time_range[2])
        t1, c1 = t1[mask1], c1[mask1]
        t2, c2 = t2[mask2], c2[mask2]
    end

    @series begin
        seriestype --> :line
        color --> colors[1]
        label --> "$(labels[1]) ($(length(t1)) pts)"
        t1, c1
    end

    @series begin
        seriestype --> :line
        color --> colors[2]
        label --> "$(labels[2]) ($(length(t2)) pts)"
        t2, c2
    end

    return Float64[], Float64[]
end
"""
# Usage
```julia
plot(cs, Val(:noise_timeline))
```

# Notes
- Shows power and S/N vs frequency in 2x1 layout
- Useful for noise property visualization
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:noise_timeline}) where {T}
    props = noise_properties(cs)

    # Create a 2x1 layout
    layout --> (2, 1)

    # Subplot 1: Power vs Frequency
    @series begin
        subplot := 1
        title := "Power Spectrum"
        xlabel := "Frequency"
        ylabel := "Power"
        seriestype := :scatter
        markersize := 2
        label := "Power"
        cs.freq, abs.(cs.power)
    end

    # Subplot 2: S/N vs Frequency
    @series begin
        subplot := 2
        title := "Signal-to-Noise Ratio"
        xlabel := "Frequency"
        ylabel := "S/N"
        seriestype := :scatter
        markersize := 2
        label := "S/N"
        cs.freq, signal_to_noise_ratio(cs)
    end
end
"""
# Usage
```julia
plot(cs, Val(:noise_properties))
```

# Notes
- Bar chart showing key noise metrics
- Displays mean S/N, noise/signal ratio, and significant detections
"""
@recipe function f(cs::CrossSpectrum{T}, ::Val{:noise_properties}) where {T}

    title --> "Noise Analysis Summary"
    xlabel --> "Property"
    ylabel --> "Value"
    grid --> true

    props = noise_properties(cs)

    # Select key properties to plot
    prop_names = ["mean_snr", "noise_to_signal_ratio", "significant_detections"]
    prop_values = [props[name] for name in prop_names]
    prop_labels = ["Mean S/N", "Noise/Signal", "Significant Bins"]

    seriestype --> :bar
    color --> [:blue, :red, :green]
    label --> ""

    return 1:length(prop_values), prop_values
end