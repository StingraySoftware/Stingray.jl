#wrapper to avoid recpies conflicts
struct BTIAnalysisPlot{T}
    eventlist::EventList{T}
end

"""
    btianalysis(events::EventList{T}) -> BTIAnalysisPlot{T}

Create a BTIAnalysisPlot object that can be plotted to show Bad Time Interval (BTI) length distribution.

This function creates a plottable object that analyzes the distribution of bad time intervals by extracting 
Good Time Intervals (GTIs) from the event metadata and computing the complementary BTIs. When plotted, it 
generates a logarithmic histogram showing the frequency distribution of BTI durations.

# Arguments
- `events::EventList{T}`: The input event list containing timing data and GTI metadata

# Returns
- `BTIAnalysisPlot{T}`: A plottable object containing the event list data

# Usage
The returned object can be plotted using the standard `plot()` function with various customization options:

```julia
# Basic BTI analysis plot
bti_plot = btianalysis(events)
plot(bti_plot, bti_analysis=true)

# Custom binning and axis limits  
plot(bti_plot, bti_analysis=true, nbins=50, xlims_range=(1e-4, 1e4))

# Analysis with custom y-axis range
plot(bti_plot, bti_analysis=true, ylims_range=(1, 1000))
```

---

    plot(bti_plot::BTIAnalysisPlot{T}; bti_analysis=false, nbins=30, min_length=1e-3, max_length=10000.0, xlims_range=nothing, ylims_range=nothing) where T

Plot a histogram of Bad Time Interval (BTI) lengths from a BTIAnalysisPlot object.

# Plot Arguments
- `bti_plot::BTIAnalysisPlot{T}`: The BTI analysis object created by `btianalysis()`
- `bti_analysis::Bool=false`: Enable BTI analysis (plot returns `nothing` if `false`)
- `nbins::Int=30`: Number of histogram bins for the distribution plot
- `min_length::Float64=1e-3`: Minimum BTI length threshold (currently unused in filtering)
- `max_length::Float64=10000.0`: Maximum BTI length threshold (currently unused in filtering)  
- `xlims_range=nothing`: Custom x-axis limits as a tuple `(min, max)`, or `nothing` for auto-scaling
- `ylims_range=nothing`: Custom y-axis limits as a tuple `(min, max)`, or `nothing` for auto-scaling

# Returns
- `Vector{Float64}`: Array of BTI lengths in seconds when BTIs are found
- `Vector{Float64}`: Dummy data `[0.5]` when no BTIs exist (creates informational plot)
- `nothing`: When `bti_analysis=false`

# Behavior
- Returns `nothing` immediately if `bti_analysis=false`
- Creates informational plots with status messages when no GTIs or BTIs are found
- Generates logarithmic histogram of BTI length distribution when BTIs exist
- Prints diagnostic statistics including total exposure and BTI lengths
- Handles various error conditions gracefully with fallback plots

# Plot Properties
- Uses logarithmic scaling on both axes for better visualization of wide dynamic ranges
- Automatic bin range calculation based on data extent
- Steel blue fill color with transparency
- Grid enabled for easier reading

# Complete Examples
```julia
using Plots

# Read event data
events = readevents("your_file.fits")

# Create BTI analysis object and plot
bti_plot = btianalysis(events)
plot(bti_plot, bti_analysis=true)

# Customized analysis with specific parameters
plot(bti_plot, 
     bti_analysis=true, 
     nbins=50, 
     xlims_range=(1e-4, 1e4),
     ylims_range=(1, 1000))

# Save the plot
savefig("bti_analysis.png")
```

# Workflow
1. Create EventList with GTI metadata: `events = readevents("file")`
2. Create BTI analysis object: `bti_plot = btianalysis(events)`
3. Generate plot: `plot(bti_plot, bti_analysis=true)`

# Notes
- Requires GTI information in `events.meta.gti` to function properly
- Short BTIs (< 1.0 second) are tracked separately in diagnostic output
- Function includes extensive error handling for missing or invalid GTI data
- The `bti_analysis=true` parameter must be set to generate the actual analysis plot
"""

# Constructor function (typically defined elsewhere)
btianalysis(events::EventList{T}) where T = BTIAnalysisPlot(events)

@recipe function f(events::BTIAnalysisPlot{T}; 
                  bti_analysis=false,
                  nbins=30,
                  min_length=1e-3,
                  max_length=10000.0,
                  xlims_range=nothing,
                  ylims_range=nothing) where T
    
    !bti_analysis && return nothing
    
    gti_event = events.eventlist
    
    function create_no_bti_plot()
        title --> "Bad Time Interval Analysis: No BTIs Found"
        xlabel --> "Observation Quality"
        ylabel --> "Coverage Status"
        grid --> false
        legend --> false
        size --> (600, 200)
        seriestype --> :scatter
        markersize --> 0
        xlims --> (0, 1)
        ylims --> (0, 1)
        annotations --> [(0.5, 0.5, "No Bad Time Intervals Found\nAll observation time covered by GTIs")]
        return [0.5], [0.5]
    end
    
    gtis = gti_event.meta.gti
    if isnothing(gtis) || isempty(gtis)
        error("No GTI information available in EventList")
    end
    
    total_exposure = gti_exposure(gti_event)
    
    if isempty(gtis) || size(gtis, 1) == 0
        println("No BTI found - all GTI")
        return create_no_bti_plot()
    end
    
    start_time = isempty(gti_event.times) ? 0.0 : minimum(gti_event.times)
    stop_time = isempty(gti_event.times) ? 0.0 : maximum(gti_event.times)
    
    btis = nothing
    try
        btis = get_btis(gtis, start_time, stop_time)
    catch
        println("No BTI found - all GTI")
        return create_no_bti_plot()
    end
    
    if isnothing(btis) || isempty(btis) || size(btis, 1) == 0
        println("No BTI found - all GTI")
        return create_no_bti_plot()
    end
    
    bti_lengths = nothing
    try
        bti_lengths = get_gti_lengths(btis)
    catch
        println("No BTI found - all GTI")
        return create_no_bti_plot()
    end
    
    if isnothing(bti_lengths) || isempty(bti_lengths) || length(bti_lengths) == 0
        println("No BTI found - all GTI")
        return create_no_bti_plot()
    end
    
    total_bti_length = 0.0
    try
        total_bti_length = get_total_gti_length(btis)
    catch
        total_bti_length = 0.0
    end
    
    # Calculate short BTI length (< 1.0 second)
    total_short_bti_length = 0.0
    try
        short_bti_mask = bti_lengths .< 1.0
        if any(short_bti_mask)
            short_btis = btis[short_bti_mask, :]
            total_short_bti_length = get_total_gti_length(short_btis)
        end
    catch
        total_short_bti_length = 0.0
    end
    
    # Print diagnostic statistics
    println("Total exposure: $(total_exposure)")
    println("Total BTI length: $(total_bti_length)")
    println("Total BTI length (short BTIs): $(total_short_bti_length)")
    
    data_min, data_max = 0.0, 1.0
    try
        data_min = minimum(bti_lengths)
        data_max = maximum(bti_lengths)
    catch
        println("Error calculating data range - no BTI found")
        return create_no_bti_plot()
    end
    
    bin_min = max(1e-6, min(1e-3, data_min * 0.9))  # Never go below 1e-6
    bin_max = max(1e-3, max(10000.0, data_max * 1.1))  # Ensure bin_max > bin_min

    # Additional safety check
    if bin_min <= 0 || !isfinite(bin_min)
        bin_min = 1e-6
    end
    if bin_max <= bin_min || !isfinite(bin_max)
        bin_max = max(1e-3, bin_min * 1000)
    end

    title --> "Distribution of Bad Time Interval Lengths"
    xlabel --> "Length of bad time interval (seconds)"
    ylabel --> "Number of intervals"
    xscale --> :log10
    yscale --> :log10
    
    if !isnothing(xlims_range)
        xlims --> xlims_range
    else
        xlims --> (bin_min, bin_max)
    end
    
    if !isnothing(ylims_range)
        ylims --> ylims_range
    else
        ylims --> (0.5, max(100, length(bti_lengths)))
    end 
    
    grid --> true
    legend --> false
    size --> (600, 400)
    
    # Manual binning for clean bars on log scale
    log_min = log10(bin_min)
    log_max = log10(bin_max)
    bin_edges = 10 .^ range(log_min, log_max, length=nbins+1)
    
    # Calculate histogram counts
    hist_counts = zeros(Int, nbins)
    for val in bti_lengths
        bin_idx = searchsortedlast(bin_edges, val)
        if bin_idx > 0 && bin_idx <= nbins
            hist_counts[bin_idx] += 1
        end
    end
    
    # Calculate bin centers (geometric mean for log scale)
    bin_centers = sqrt.(bin_edges[1:end-1] .* bin_edges[2:end])
    
    # Filter out empty bins
    valid_bins = hist_counts .> 0
    plot_x = bin_centers[valid_bins]
    plot_y = hist_counts[valid_bins]
    
    seriestype --> :bar
    fillcolor --> :steelblue
    fillalpha --> 0.8
    linecolor --> :steelblue
    linewidth --> 1
    
    return plot_x, plot_y
end
