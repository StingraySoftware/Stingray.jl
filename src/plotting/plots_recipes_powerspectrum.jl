"""
    @recipe f(ps::AbstractPowerSpectrum{T}; kwargs...)

Recipe for plotting power spectra with comprehensive customization options.

# Keyword Arguments
- `show_noise=false`: Display Poisson noise level
- `show_errors=false`: Show error bars if available
- `freq_mult=false`: Multiply power by frequency (ν×P vs P)
- `log_scale=true`: Use logarithmic axes
- `noise_style=:dash`: Line style for noise level
- `meanrate=nothing`: Mean count rate (uses stored value if available)
- `segment_size=nothing`: Segment size information
- `drawstyle=:steppost`: Plot drawing style
- `subtract_noise=false`: Remove noise level from power
- `error_alpha=0.3`: Transparency for error bars

# Examples
```julia
ps = AveragedPowerspectrum(events, 256.0, norm="frac")
plot(ps, show_noise=true, freq_mult=true)
```
"""
@recipe function f(ps::AbstractPowerSpectrum{T};
                  show_noise=false,
                  show_errors=false,
                  freq_mult=false,
                  log_scale=true,
                  noise_style=:dash,
                  meanrate=nothing,
                  segment_size=nothing,
                  drawstyle=:steppost,
                  subtract_noise=false,
                  error_alpha=0.3) where T
    
    isempty(ps.freq) && error("PowerSpectrum is empty")
    
    freq_vec = ps.freq
    power_vec = copy(ps.power)
    
    effective_meanrate = if !isnothing(meanrate)
        meanrate
    elseif hasfield(typeof(ps), :mean_rate)
        ps.mean_rate
    else
        nothing
    end
    
    if subtract_noise && show_noise && !isnothing(effective_meanrate)
        noise_level = poisson_level(ps.norm; meanrate=effective_meanrate)
        power_vec .-= noise_level
        if log_scale
            power_vec = max.(power_vec, 1e-10 * maximum(power_vec))
        end
    end
    
    if freq_mult
        power_vec .*= freq_vec
        ylabel_text = get_norm_label(ps.norm) * " × f"
    else
        ylabel_text = get_norm_label(ps.norm)
    end
    
    noise_level = nothing
    if show_noise && !isnothing(effective_meanrate)
        noise_level = get_poisson_level(ps.norm; meanrate=effective_meanrate)
        if freq_mult
            noise_level *= mean(freq_vec)
        end
    end
    
    title --> "Power Spectrum"
    xlabel --> "Frequency (Hz)"
    ylabel --> ylabel_text
    grid --> true
    minorgrid --> true  
    legend --> :topright
    
    if log_scale
        xscale --> :log10
        yscale --> :log10
    end
    
    if show_errors && !isnothing(ps.power_err)
        error_vec = copy(ps.power_err)
        if freq_mult
            error_vec .*= freq_vec
        end
        
        @series begin
            seriestype := :scatter
            marker := :none
            yerror := error_vec
            errorbar_color := :gray
            errorbar_width := 0.5
            errorbar_alpha := error_alpha
            label := ""
            freq_vec, power_vec
        end
    end
    
    if show_noise && !isnothing(noise_level)
        @series begin
            seriestype := :hline
            linestyle := noise_style
            linewidth := 2
            color := :red
            label := subtract_noise ? "Original Noise Level" : "Poisson Noise Level"
            y := [noise_level]
        end
    end
    
    seriestype --> drawstyle
    linewidth --> 1.5
    color --> :blue
    label --> subtract_noise ? "Noise-Subtracted Power" : "Power Spectrum"
    
    return freq_vec, power_vec
end

"""
    @recipe f(events::EventList{Vector{T}, M}; kwargs...)

Recipe for plotting power spectra directly from event lists.

# Keyword Arguments
- `segment_size=256.0`: Length of segments for averaging (seconds)
- `dt=0.001`: Time resolution for events (seconds)
- `norm="frac"`: Normalization type
- `show_noise=true`: Display Poisson noise level
- `show_errors=false`: Show error bars
- `freq_mult=false`: Multiply power by frequency
- `log_scale=true`: Use logarithmic axes
- `drawstyle=:steppost`: Plot drawing style

# Examples
```julia
events = readevents("data.evt")
plot(events, segment_size=256.0, norm="frac", show_noise=true)
```
"""
@recipe function f(events::EventList{Vector{T}, M};
                  segment_size=256.0,
                  dt=0.001,
                  norm="frac",
                  show_noise=true,
                  show_errors=false,
                  freq_mult=false,
                  log_scale=true,
                  drawstyle=:steppost) where {T<:Real, M}
    
    isempty(events.times) && error("EventList is empty")
    
    ps = AveragedPowerspectrum(events, segment_size; norm=norm, dt=dt)
    
    @series begin
        show_noise := show_noise
        show_errors := show_errors
        freq_mult := freq_mult
        log_scale := log_scale
        drawstyle := drawstyle
        ps
    end
end

"""
    @recipe f(lc::LightCurve{T}; kwargs...)

Recipe for plotting power spectra from light curves.

# Keyword Arguments
- `segment_size=nothing`: Segment size (auto-determined if nothing)
- `norm="frac"`: Normalization type
- `show_noise=true`: Display Poisson noise level
- `show_errors=false`: Show error bars
- `freq_mult=false`: Multiply power by frequency
- `log_scale=true`: Use logarithmic axes
- `drawstyle=:steppost`: Plot drawing style

# Examples
```julia
lc = create_lightcurve(events, 0.1)
plot(lc, segment_size=128, norm="frac", show_noise=true)
```
"""
@recipe function f(lc::LightCurve{T};
                  segment_size=nothing,
                  norm="frac",
                  show_noise=true,
                  show_errors=false,
                  freq_mult=false,
                  log_scale=true,
                  drawstyle=:steppost) where T
    
    isempty(lc.time) && error("LightCurve is empty")
    
    if isnothing(segment_size)
        total_time = lc.metadata.time_range[2] - lc.metadata.time_range[1]
        segment_size = total_time / 32
    end
    
    ps = if segment_size < (length(lc.time) * lc.metadata.bin_size) / 2
        AveragedPowerspectrum(lc, segment_size; norm=norm)
    else
        Powerspectrum(lc; norm=norm) 
    end
    
    power_vec = copy(ps.power)
    
    if freq_mult
        power_vec .*= ps.freq
        ylabel_text = get_norm_label(ps.norm) * " × f"
    else
        ylabel_text = get_norm_label(ps.norm)
    end
    
    title --> "Power Spectrum from Light Curve"
    xlabel --> "Frequency (Hz)"
    ylabel --> ylabel_text
    grid --> true
    minorgrid --> true  
    legend --> :topright
    
    if log_scale
        xscale --> :log10
        yscale --> :log10
    end
    
    if show_errors && !isnothing(ps.power_err)
        error_vec = copy(ps.power_err)
        if freq_mult
            error_vec .*= ps.freq
        end
        
        @series begin
            seriestype := :scatter
            marker := :none
            yerror := error_vec
            errorbar_color := :gray  
            errorbar_width := 0.5
            errorbar_alpha := 0.3
            label := ""
            ps.freq, power_vec
        end
    end
    
    if show_noise
        effective_meanrate = if hasfield(typeof(ps), :mean_rate)
            ps.mean_rate
        else
            total_counts = sum(lc.counts)
            total_time = lc.metadata.time_range[2] - lc.metadata.time_range[1]
            total_counts / total_time
        end
        
        noise_level = get_poisson_level(ps.norm; meanrate=effective_meanrate)
        if freq_mult
            noise_level *= mean(ps.freq)
        end
        
        @series begin
            seriestype := :hline
            linestyle := :dash
            linewidth := 2
            color := :red
            label := "Poisson Noise Level"
            y := [noise_level]
        end
    end
    
    seriestype --> drawstyle
    linewidth --> 1.5
    color --> :blue
    label --> "Power Spectrum"
    
    return ps.freq, power_vec
end

"""
    @recipe f(spectra::Vector{<:AbstractPowerSpectrum}; kwargs...)

Recipe for comparing multiple power spectra on the same plot.

# Keyword Arguments
- `labels=nothing`: Custom labels for each spectrum
- `colors=nothing`: Custom colors for each spectrum
- `show_noise=false`: Display noise levels
- `freq_mult=false`: Multiply power by frequency
- `log_scale=true`: Use logarithmic axes
- `drawstyle=:steppost`: Plot drawing style
- `alpha=1.0`: Line transparency
- `linewidth=1.5`: Line thickness

# Examples
```julia
ps1 = AveragedPowerspectrum(events, 128.0, norm="frac")
ps2 = AveragedPowerspectrum(events, 256.0, norm="frac")
plot([ps1, ps2], labels=["128s", "256s"])
```
"""
@recipe function f(spectra::Vector{<:AbstractPowerSpectrum};
                  labels=nothing,
                  colors=nothing,
                  show_noise=false,
                  freq_mult=false,
                  log_scale=true,
                  drawstyle=:steppost,
                  alpha=1.0,
                  linewidth=1.5)
    
    isempty(spectra) && error("Spectra vector is empty")
    
    title --> "Power Spectra Comparison"
    xlabel --> "Frequency (Hz)"
    ylabel --> freq_mult ? get_norm_label(spectra[1].norm) * " × f" : get_norm_label(spectra[1].norm)
    grid --> true
    minorgrid --> true
    legend --> :topright
    
    if log_scale
        xscale --> :log10
        yscale --> :log10
    end
    
    default_colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :cyan, :magenta]
    plot_colors = isnothing(colors) ? default_colors : colors
    
    for (i, ps) in enumerate(spectra)
        power_vec = freq_mult ? ps.power .* ps.freq : ps.power
        
        @series begin
            seriestype := drawstyle
            linewidth := linewidth
            color := plot_colors[mod1(i, length(plot_colors))]
            alpha := alpha
            label := isnothing(labels) ? "Spectrum $i" : labels[i]
            
            ps.freq, power_vec
        end
    end
end

"""
    @recipe f(events_dict::Dict{String, EventList}; kwargs...)

Recipe for multi-band spectral-timing analysis.

# Keyword Arguments
- `segment_size=256.0`: Segment size for averaging
- `norm="frac"`: Normalization type
- `energy_labels=nothing`: Custom labels for energy bands
- `colors=nothing`: Custom colors for each band
- `freq_mult=false`: Multiply power by frequency
- `log_scale=true`: Use logarithmic axes

# Examples
```julia
events_dict = Dict(
    "0.5-2 keV" => events_soft,
    "2-10 keV" => events_hard
)
plot(events_dict, segment_size=256.0, energy_labels=["Soft", "Hard"])
```
"""
@recipe function f(events_dict::Dict{String, EventList};
                  segment_size=256.0,
                  norm="frac",
                  energy_labels=nothing,
                  colors=nothing,
                  freq_mult=false,
                  log_scale=true)
    
    title --> "Multi-Band Power Spectra"
    xlabel --> "Frequency (Hz)"
    ylabel --> freq_mult ? get_norm_label(norm) * " × f" : get_norm_label(norm)
    
    if log_scale
        xscale --> :log10
        yscale --> :log10
    end
    
    default_colors = [:blue, :red, :green, :orange, :purple, :brown]
    plot_colors = isnothing(colors) ? default_colors : colors
    
    band_names = collect(keys(events_dict))
    labels = isnothing(energy_labels) ? band_names : energy_labels
    
    for (i, (band_name, events)) in enumerate(events_dict)
        ps = AveragedPowerspectrum(events, segment_size; norm=norm)
        power_vec = freq_mult ? ps.power .* ps.freq : ps.power
        
        @series begin
            seriestype := :steppost
            linewidth := 1.5
            color := plot_colors[mod1(i, length(plot_colors))]
            label := labels[i]
            
            ps.freq, power_vec
        end
    end
end