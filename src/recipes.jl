"""
    f(el::EventList{T}, bin_size::Real=0.1; tstart=nothing, tseg=nothing, show_errors=true, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T

Recipe for plotting an EventList directly as a light curve.

# Arguments
- `el`: EventList containing time and energy data
- `bin_size`: Size of time bins for the light curve
- `tstart`: Optional starting time for light curve (defaults to minimum time)
- `tseg`: Optional time segment length (defaults to full range)
- `show_errors`: Whether to display Poisson error bars
- `show_gaps`: Whether to highlight data gaps
- `gap_threshold`: Multiplier of median time difference to identify gaps
- `axis_limits`: Optional array of [xmin, xmax, ymin, ymax] for plot limits
"""
@recipe function f(el::EventList{T}, bin_size::Real=0.1; tstart=nothing, tseg=nothing, show_errors=true, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T
    bin_size_t = convert(T, bin_size)
    
    if isempty(el.times)
        @warn "Event list contains no time data. Cannot create light curve."
        return [], []
    end
    
    min_time = minimum(el.times)
    max_time = maximum(el.times)
    
    start_time = isnothing(tstart) ? min_time : convert(T, tstart)
    end_time = isnothing(tseg) ? max_time : start_time + convert(T, tseg)
    
    filtered_times = filter(t -> start_time <= t <= end_time, el.times)
    
    edges = start_time:bin_size_t:(end_time + bin_size_t)
    bin_centers = edges[1:end-1] .+ bin_size_t/2
    n_bins = length(bin_centers)
    
    counts = zeros(Int, n_bins)
    for t in filtered_times
        bin_idx = floor(Int, (t - start_time) / bin_size_t) + 1
        if 1 <= bin_idx <= n_bins
            counts[bin_idx] += 1
        end
    end
    
    errors = sqrt.(float.(counts))
    
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits
            
            if !isnothing(xmin) && !isnothing(xmax)
                xlims --> (xmin, xmax)
            elseif !isnothing(xmin)
                xlims --> (xmin, :auto)
            elseif !isnothing(xmax)
                xlims --> (:auto, xmax)
            end
            
            if !isnothing(ymin) && !isnothing(ymax)
                ylims --> (ymin, ymax)
            elseif !isnothing(ymin)
                ylims --> (ymin, :auto)
            elseif !isnothing(ymax)
                ylims --> (:auto, ymax)
            end
        else
            @warn "axis_limits should be an array of 4 elements [xmin, xmax, ymin, ymax]"
        end
    end
    
    if show_gaps
        sorted_times = sort(filtered_times)
        time_diffs = diff(sorted_times)
        median_diff = median(time_diffs)
        threshold = gap_threshold * median_diff
        gap_indices = findall(diff -> diff > threshold, time_diffs)
        
        gap_starts = sorted_times[gap_indices]
        gap_ends = sorted_times[gap_indices .+ 1]
        
        @series begin
            title := "Light Curve with Gaps Highlighted"
            xlabel := "Time"
            ylabel := "Counts"
            seriestype := :steppost
            linewidth := 2
            color := :blue
            label := "Light Curve (bin size: $bin_size_t)"
            bin_centers, counts
        end
        
        max_count = isempty(counts) ? one(T) : maximum(counts)
        
        if !isempty(gap_starts)
            for (i, (start, stop)) in enumerate(zip(gap_starts, gap_ends))
                @series begin
                    seriestype := :shape
                    fillalpha := 0.3
                    fillcolor := :red
                    linecolor := :red
                    label := i == 1 ? "Gaps (threshold: $(round(threshold, digits=2)))" : nothing
                    [start, stop, stop, start, start], 
                    [zero(T), zero(T), max_count, max_count, zero(T)]
                end
            end
        end
        
        return [], []
    else
        title --> "Light Curve"
        xlabel --> "Time"
        ylabel --> "Counts"
        
        if show_errors
            seriestype --> :scatter
            yerror --> errors
            label --> "Light Curve with Poisson errors"
        else
            seriestype --> :steppost
            label --> "Light Curve (bin size: $bin_size_t)"
        end
        
        return bin_centers, counts
    end
end

"""
    f(lc::LightCurve{T}; show_errors=true, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T

Recipe for plotting a LightCurve object with optional error bars and gap highlighting.

# Arguments
- `lc`: LightCurve object containing binned time series data
- `show_errors`: Whether to display error bars
- `show_gaps`: Whether to highlight data gaps
- `gap_threshold`: Multiplier of median time difference to identify gaps
- `axis_limits`: Optional array of [xmin, xmax, ymin, ymax] for plot limits
"""
@recipe function f(lc::LightCurve{T}; show_errors=true, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits
            
            if !isnothing(xmin) && !isnothing(xmax)
                xlims --> (xmin, xmax)
            elseif !isnothing(xmin)
                xlims --> (xmin, :auto)
            elseif !isnothing(xmax)
                xlims --> (:auto, xmax)
            end
            
            if !isnothing(ymin) && !isnothing(ymax)
                ylims --> (ymin, ymax)
            elseif !isnothing(ymin)
                ylims --> (ymin, :auto)
            elseif !isnothing(ymax)
                ylims --> (:auto, ymax)
            end
        else
            @warn "axis_limits should be an array of 4 elements [xmin, xmax, ymin, ymax]"
        end
    end
    
    if !show_gaps
        title --> "Light Curve"
        xlabel --> "Time"
        ylabel --> "Counts"
        
        if show_errors
            seriestype --> :scatter
            yerror --> lc.count_error
            label --> "Light Curve with $(lc.err_method) errors"
        else
            seriestype --> :steppost
            label --> "Light Curve"
        end
        
        return lc.timebins, lc.counts
    else
        if length(lc.timebins) > 1
            bin_width = lc.timebins[2] - lc.timebins[1]
        else
            bin_width = one(T)
        end
        
        pseudo_events = T[]
        for (t, count) in zip(lc.timebins, lc.counts)
            if count > 0
                if count == 1
                    push!(pseudo_events, t)
                else
                    half_width = bin_width / 2
                    event_spacing = bin_width / count
                    for i in 1:count
                        event_time = t - half_width + (i-0.5)*event_spacing
                        push!(pseudo_events, event_time)
                    end
                end
            end
        end
        
        sort!(pseudo_events)
        
        gaps = Tuple{T,T}[]
        if length(pseudo_events) > 1
            time_diffs = diff(pseudo_events)
            median_diff = median(time_diffs)
            threshold = gap_threshold * median_diff
            
            gap_indices = findall(diff -> diff > threshold, time_diffs)
            
            for idx in gap_indices
                start_time = pseudo_events[idx]
                end_time = pseudo_events[idx+1]
                push!(gaps, (start_time, end_time))
            end
        end
        
        @series begin
            title := "Light Curve with Gaps Highlighted"
            xlabel := "Time"
            ylabel := "Counts"
            seriestype := :steppost
            linewidth := 2
            color := :blue
            label := "Light Curve"
            lc.timebins, lc.counts
        end
        
        max_count = isempty(lc.counts) ? one(T) : maximum(lc.counts)
        
        if !isempty(gaps)
            for (i, (start, stop)) in enumerate(gaps)
                @series begin
                    seriestype := :shape
                    fillalpha := 0.3
                    fillcolor := :red
                    linecolor := :red
                    label := i == 1 ? "Gaps" : nothing
                    [start, stop, stop, start, start], 
                    [zero(T), zero(T), max_count, max_count, zero(T)]
                end
            end
        end
        
        return [], []
    end
end

"""
    f(::Type{LightCurve}, el::EventList{T}, bin_size::Real=0.1; 
      tstart=nothing, tseg=nothing, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T

Recipe that converts an EventList to a LightCurve for plotting.

# Arguments
- `el`: EventList containing time and energy data
- `bin_size`: Size of time bins for the light curve
- `tstart`: Optional starting time for light curve
- `tseg`: Optional time segment length
- `show_gaps`: Whether to highlight data gaps
- `gap_threshold`: Multiplier of median time difference to identify gaps
- `axis_limits`: Optional array of [xmin, xmax, ymin, ymax] for plot limits
"""
@recipe function f(::Type{LightCurve}, el::EventList{T}, bin_size::Real=0.1; 
                  tstart=nothing, tseg=nothing, show_gaps=false, gap_threshold=10.0, axis_limits=nothing) where T
    bin_size_t = convert(T, bin_size)
    
    min_time = minimum(el.times)
    max_time = maximum(el.times)
    
    start_time = isnothing(tstart) ? min_time : convert(T, tstart)
    end_time = isnothing(tseg) ? max_time : start_time + convert(T, tseg)
    
    filtered_times = filter(t -> start_time <= t <= end_time, el.times)
    
    edges = start_time:bin_size_t:(end_time + bin_size_t)
    bin_centers = edges[1:end-1] .+ bin_size_t/2
    n_bins = length(bin_centers)
    
    counts = zeros(Int, n_bins)
    for t in filtered_times
        bin_idx = floor(Int, (t - start_time) / bin_size_t) + 1
        if 1 <= bin_idx <= n_bins
            counts[bin_idx] += 1
        end
    end
    
    errors = sqrt.(float.(counts))
    
    lc = LightCurve{T}(
        bin_centers,
        convert.(T, counts),
        convert.(T, errors),
        :poisson
    )
    
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits
            
            if !isnothing(xmin) && !isnothing(xmax)
                xlims --> (xmin, xmax)
            elseif !isnothing(xmin)
                xlims --> (xmin, :auto)
            elseif !isnothing(xmax)
                xlims --> (:auto, xmax)
            end
            
            if !isnothing(ymin) && !isnothing(ymax)
                ylims --> (ymin, ymax)
            elseif !isnothing(ymin)
                ylims --> (ymin, :auto)
            elseif !isnothing(ymax)
                ylims --> (:auto, ymax)
            end
        else
            @warn "axis_limits should be an array of 4 elements [xmin, xmax, ymin, ymax]"
        end
    end
    
    if show_gaps
        sorted_times = sort(filtered_times)
        time_diffs = diff(sorted_times)
        median_diff = median(time_diffs)
        threshold = gap_threshold * median_diff
        gap_indices = findall(diff -> diff > threshold, time_diffs)
        
        gap_starts = sorted_times[gap_indices]
        gap_ends = sorted_times[gap_indices .+ 1]
        
        @series begin
            title := "Light Curve with Gaps Highlighted"
            xlabel := "Time"
            ylabel := "Counts"
            seriestype := :steppost
            linewidth := 2
            color := :blue
            label := "Light Curve (bin size: $bin_size_t)"
            lc.timebins, lc.counts
        end
        
        max_count = isempty(counts) ? one(T) : maximum(counts)
        
        if !isempty(gap_starts)
            for (i, (start, stop)) in enumerate(zip(gap_starts, gap_ends))
                @series begin
                    seriestype := :shape
                    fillalpha := 0.3
                    fillcolor := :red
                    linecolor := :red
                    label := i == 1 ? "Gaps (threshold: $(round(threshold, digits=2)))" : nothing
                    [start, stop, stop, start, start], 
                    [zero(T), zero(T), max_count, max_count, zero(T)]
                end
            end
        end
        
        return [], []
    else
        title --> "Light Curve"
        xlabel --> "Time"
        ylabel --> "Counts"
        seriestype --> :steppost
        label --> "Light Curve (bin size: $bin_size_t)"
        
        return lc.timebins, lc.counts
    end
end

"""
    f(::Type{Val{:events}}, x, y, z; color=:auto, markersize=4)

Recipe for rendering a general events timeline using vertical markers.
Creates a scatter plot with vertical line markers positioned at the x-coordinates.

# Arguments
- `color`: Color for the event markers
- `markersize`: Size of the event markers
"""
@recipe function f(::Type{Val{:events}}, x, y, z; color=:auto, markersize=4)
    seriestype := :scatter
    markershape := :vline
    markersize := markersize
    markerstrokecolor := color
    label := "Events"
    
    y = ones(length(x))
    
    x, y
end

"""
    f(el::EventList{T}, ::Type{Val{:events}}; 
      tstart=nothing, tseg=nothing, axis_limits=nothing, color_by_energy=false) where T

Recipe for rendering an EventList as a timeline visualization.

# Arguments
- `tstart`: Optional starting time for the plot window
- `tseg`: Optional time segment length
- `axis_limits`: Optional array of [xmin, xmax, ymin, ymax] for plot limits
- `color_by_energy`: If true, colors the events by their energy values
"""
@recipe function f(el::EventList{T}, ::Type{Val{:events}}; 
                  tstart=nothing, tseg=nothing, axis_limits=nothing, color_by_energy=false) where T
    min_time = minimum(el.times)
    max_time = maximum(el.times)
    
    start_time = isnothing(tstart) ? min_time : convert(T, tstart)
    end_time = isnothing(tseg) ? max_time : start_time + convert(T, tseg)
    
    time_indices = findall(t -> start_time <= t <= end_time, el.times)
    filtered_times = el.times[time_indices]
    
    if !isnothing(axis_limits) && length(axis_limits) == 4
        xmin, xmax, ymin, ymax = axis_limits
        
        if !isnothing(xmin) && !isnothing(xmax)
            xlims --> (xmin, xmax)
        elseif !isnothing(xmin)
            xlims --> (xmin, :auto)
        elseif !isnothing(xmax)
            xlims --> (:auto, xmax)
        end
        
        if !isnothing(ymin) && !isnothing(ymax)
            ylims --> (ymin, ymax)
        elseif !isnothing(ymin)
            ylims --> (ymin, :auto)
        elseif !isnothing(ymax)
            ylims --> (:auto, ymax)
        end
    end
    
    title --> "Event Timeline"
    xlabel --> "Time"
    ylabel --> ""
    yticks --> []
    grid --> false
    
    if color_by_energy
        filtered_energies = el.energies[time_indices]
        
        seriestype := :scatter
        markershape := :vline
        markersize := 10
        marker_z := filtered_energies
        colorbar --> true
        label := "Events (colored by energy)"
        
        return filtered_times, zeros(length(filtered_times))
    else
        seriestype := :events
        return filtered_times, []
    end
end

"""
    f(::Type{Val{:lightcurve_with_events}}, lc::LightCurve{T}, el::EventList{T}; 
      bin_size=0.1, tstart=nothing, tseg=nothing, axis_limits=nothing, color_by_energy=false) where T

Recipe for rendering a combined visualization of a LightCurve with events from an EventList.

# Arguments
- `bin_size`: Bin size for light curve (if generating from events)
- `tstart`: Optional starting time for the plot window
- `tseg`: Optional time segment length
- `axis_limits`: Optional array of [xmin, xmax, ymin, ymax] for plot limits
- `color_by_energy`: If true, colors the events by their energy values
"""
@recipe function f(::Type{Val{:lightcurve_with_events}}, lc::LightCurve{T}, el::EventList{T}; 
                  bin_size=0.1, tstart=nothing, tseg=nothing, axis_limits=nothing, 
                  color_by_energy=false) where T
    
    @series begin
        seriestype := :steppost
        label := "Light Curve"
        lc.timebins, lc.counts
    end
    
    min_time = isnothing(tstart) ? minimum(el.times) : convert(T, tstart)
    max_time = isnothing(tseg) ? maximum(el.times) : min_time + convert(T, tseg)
    
    time_indices = findall(t -> min_time <= t <= max_time, el.times)
    filtered_times = el.times[time_indices]
    
    @series begin
        seriestype := :scatter
        markershape := :vline
        markersize := 10
        
        if color_by_energy
            filtered_energies = el.energies[time_indices]
            marker_z := filtered_energies
            colorbar := true
            label := "Events (energy)"
        else
            markerstrokecolor := :red
            label := "Events"
        end
        
        filtered_times, zeros(length(filtered_times))
    end
    
    if !isnothing(axis_limits) && length(axis_limits) == 4
        xmin, xmax, ymin, ymax = axis_limits
        
        if !isnothing(xmin) && !isnothing(xmax)
            xlims --> (xmin, xmax)
        elseif !isnothing(xmin)
            xlims --> (xmin, :auto)
        elseif !isnothing(xmax)
            xlims --> (:auto, xmax)
        end
        
        if !isnothing(ymin) && !isnothing(ymax)
            ylims --> (ymin, ymax)
        elseif !isnothing(ymin)
            ylims --> (ymin, :auto)
        elseif !isnothing(ymax)
            ylims --> (:auto, ymax)
        end
    end
    
    title --> "Light Curve with Events"
    xlabel --> "Time"
    ylabel --> "Counts"
    
    return [], []
end