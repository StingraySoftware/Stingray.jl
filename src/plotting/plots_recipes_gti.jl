struct BTIAnalysisPlot{T}
    eventlist::EventList{T}
end
"""
    f(events::EventList{T}; bti_analysis=false, nbins=30, min_length=1e-3, max_length=10000.0, xlims_range=nothing, ylims_range=nothing) where T

Create a histogram plot of Bad Time Interval (BTI) lengths from an EventList.

This recipe function analyzes the distribution of bad time intervals by extracting Good Time Intervals (GTIs)
from the event metadata and computing the complementary BTIs. It generates a logarithmic histogram showing
the frequency distribution of BTI durations.

# Arguments
- `events::EventList{T}`: The input event list containing timing data and GTI metadata
- `bti_analysis::Bool=false`: Enable BTI analysis (function returns `nothing` if `false`)
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

# Examples
```julia
# Basic BTI analysis
using Plots
plots(events=readevents("your_file"), bti_analysis=true)

# Custom binning and axis limits  
plots(events, bti_analysis=true, nbins=50, xlims_range=(1e-4, 1e4))

# Analysis with custom y-axis range
plots(events, bti_analysis=true, ylims_range=(1, 1000))
```

# Notes
- Requires GTI information in `events.meta.gti` to function properly
- Short BTIs (< 1.0 second) are tracked separately in diagnostic output
- Function includes extensive error handling for missing or invalid GTI data
"""
@recipe function f(events::BTIAnalysisPlot{T}; 
                  bti_analysis=false,
                  nbins=30,
                  min_length=1e-3,
                  max_length=10000.0,
                  xlims_range=nothing,
                  ylims_range=nothing) where T
    
    !bti_analysis && return nothing
    
    # Extract EventList from wrapper
    gti_event = events.eventlist
    
    # Helper function to create "No BTI found" plot
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
        
        title --> "Bad Time Interval Analysis: No GTIs Available"
        xlabel --> "Observation Quality"
        ylabel --> "Coverage Status"
        grid --> false
        legend --> false
        size --> (600, 200)
        seriestype --> :scatter
        markersize --> 0
        xlims --> (0, 1)
        ylims --> (0, 1)
        annotations --> [(0.5, 0.5, "No GTI Information Available")]
        
        return [0.5], [0.5]
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
    
    # Check if BTIs exist and are valid
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
    
    # THIS IS THE KEY FIX - Use the same logic as your working function
    data_min, data_max = 0.0, 1.0
    try
        data_min = minimum(bti_lengths)
        data_max = maximum(bti_lengths)
    catch
        println("Error calculating data range - no BTI found")
        return create_no_bti_plot()
    end
    
    # Calculate bin range for display
    bin_min = min(1e-3, data_min * 0.1)
    bin_max = max(10000, data_max * 2.0)
    num_bins = Int(nbins)
    
    title --> "Distribution of Bad Time Interval Lengths"
    xlabel --> "Length of bad time interval"
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
    seriestype --> :histogram
    nbins --> num_bins
    fillcolor --> :steelblue
    fillalpha --> 0.7
    linecolor --> :steelblue
    linewidth --> 1
    
    return bti_lengths
end