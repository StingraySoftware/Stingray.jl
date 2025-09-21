using Plots:text
function create_segments(events::EventList, segment_duration::Real; bin_size::Real = 1.0)
    # Get actual time range from the data
    t_start = minimum(events.times)
    t_stop = maximum(events.times)
    total_time = t_stop - t_start
    println("Data time range: $(t_start) to $(t_stop) ($(total_time) seconds total)")

    # Calculate number of segments
    n_segments = ceil(Int, total_time / segment_duration)
    println("Creating $(n_segments) segments of $(segment_duration) seconds each")

    segments = Vector{LightCurve}()

    for i = 1:n_segments
        start_time = t_start + (i - 1) * segment_duration
        stop_time = min(t_start + i * segment_duration, t_stop)

        # Filter events in this segment using matrix operations
        event_times = events.times
        mask = (event_times .>= start_time) .& (event_times .<= stop_time)
        segment_events = event_times[mask]

        events_in_segment = sum(mask)
        println(
            "Segment $i: $(events_in_segment) events from $(start_time) to $(stop_time)",
        )

        # Create time bins for this segment
        time_edges = collect(start_time:bin_size:stop_time)
        if length(time_edges) < 2
            time_edges = [start_time, stop_time]
        end

        # Bin centers
        time_centers = (time_edges[1:end-1] + time_edges[2:end]) / 2
        n_bins = length(time_centers)

        # Initialize counts
        counts = zeros(Int, n_bins)

        # Histogram events into bins using matrix operations
        if events_in_segment > 0
            # Use searchsortedlast for efficient binning
            for event_time in segment_events
                bin_idx = searchsortedlast(time_edges, event_time)
                if bin_idx > 0 && bin_idx <= n_bins
                    counts[bin_idx] += 1
                end
            end
        end

        # Create light curve data matrix: [time, counts, errors]
        # Calculate Poisson errors
        count_errors = sqrt.(max.(counts, 1))  # Avoid sqrt(0)

        # Create the data matrix
        lc_matrix = hcat(time_centers, counts, count_errors)

        # Create dummy metadata
        dummy_metadata = LightCurveMetadata(
            "Unknown",
            "Unknown",
            "Unknown",
            0.0,
            (start_time, stop_time),
            bin_size,
            Vector{Dict{String,Any}}(),
            Dict{String,Any}(),
        )

        # Create LightCurve object
        lc = LightCurve(
            lc_matrix[:, 1],  # time
            bin_size,         # dt
            Int.(lc_matrix[:, 2]),  # counts
            lc_matrix[:, 3],  # count_error
            nothing,          # exposure
            Vector{EventProperty{Float64}}(),  # properties
            dummy_metadata,   # metadata
            :poisson,         # err_method
        )

        push!(segments, lc)
    end

    return segments
end
"""
    plot(el::EventList{T}, bin_size::Real=1.0; kwargs...)

Plot a light curve from an EventList with optional Good Time Intervals (GTIs) and Bad Time Intervals (BTIs).

# Arguments
- `el::EventList{T}`: Event list containing photon arrival times
- `bin_size::Real=1.0`: Time bin size in seconds

# Keywords
- `tstart=nothing`: Start time for the light curve (defaults to first event)
- `tstop=nothing`: Stop time for the light curve (defaults to last event)
- `energy_filter=nothing`: Energy range filter as (min, max) tuple
- `show_errors=false`: Display error bars using specified error method
- `show_btis=false`: Show Bad Time Intervals as red shaded regions
- `show_bti=false`: Alias for `show_btis`
- `show_gtis=false`: Show Good Time Intervals as green shaded regions
- `show_gti=false`: Alias for `show_gtis`
- `gtis=nothing`: GTI matrix to use (overrides file/metadata GTIs)
- `gti_file=nothing`: FITS file containing GTI extension
- `gti_hdu="GTI"`: HDU name for GTI data
- `bti_alpha=0.3`: Transparency for BTI shading
- `gti_alpha=0.2`: Transparency for GTI shading
- `gap_threshold=10.0`: Minimum gap size to consider as BTI
- `axis_limits=nothing`: Plot limits as `[xmin, xmax]` or `[xmin, xmax, ymin, ymax]`
- `err_method=:poisson`: Error calculation method (`:poisson`, `:gaussian`)

# Returns
- `Tuple{Vector, Vector}`: Time bins and corresponding count rates

# Examples
```julia
# Basic light curve
plot(events, 1.0)

# With error bars and GTIs
plot(events, 0.5, show_errors=true, show_gtis=true)

# Custom time range with energy filter
plot(events, 2.0, tstart=100.0, tstop=500.0, energy_filter=(2.0, 10.0))

# With custom axis limits
plot(events, 1.0, axis_limits=[0, 1000, 0, 100])
```

# Notes
- GTI priority: explicit `gtis` > `gti_file` > `el.meta.gti`
- BTIs are calculated as gaps between GTIs exceeding `gap_threshold`
- Error bars use Poisson statistics by default
"""
@recipe function f(
    el::EventList{Vector{T}, M},
    bin_size::Real = 1.0;
    tstart = nothing,
    tstop = nothing,
    energy_filter = nothing,
    show_errors = false,
    show_btis = false,
    show_bti = false,
    show_gtis = false,
    show_gti = false,
    gtis = nothing,
    gti_file = nothing,
    gti_hdu = "GTI",
    bti_alpha = 0.3,
    gti_alpha = 0.2,
    gap_threshold = 10.0,
    axis_limits = nothing,
    err_method = :poisson,
) where {T<:Real, M<:FITSMetadata}

    isempty(el.times) && error("EventList is empty")
    
    show_gtis = show_gtis || show_gti
    show_btis = show_btis || show_bti

    # Apply filters directly to EventList data
    filtered_el = deepcopy(el)  # Work with a copy
    
    # Apply energy filter if specified
    if !isnothing(energy_filter) && has_energies(filtered_el)
        e_min, e_max = energy_filter
        filter_energy!(e -> e_min ≤ e ≤ e_max, filtered_el)
    end
    
    # Apply time filter if specified
    if !isnothing(tstart) || !isnothing(tstop)
        t_min = isnothing(tstart) ? minimum(times(filtered_el)) : tstart
        t_max = isnothing(tstop) ? maximum(times(filtered_el)) : tstop
        filter_time!(t -> t_min ≤ t ≤ t_max, filtered_el)
    end
    
    # Check if we have events after filtering
    if isempty(times(filtered_el))
        error("No events remain after applying filters")
    end
    
    # Get filtered times
    event_times = times(filtered_el)
    time_min = minimum(event_times)
    time_max = maximum(event_times)
    
    # Create histogram bins directly from event times
    plot_tstart = isnothing(tstart) ? time_min : max(tstart, time_min)
    plot_tstop = isnothing(tstop) ? time_max : min(tstop, time_max)
    
    # Create time bins
    n_bins = max(1, Int(ceil((plot_tstop - plot_tstart) / bin_size)))
    bin_edges = range(plot_tstart, plot_tstop, length=n_bins+1)
    bin_centers = (bin_edges[1:end-1] + bin_edges[2:end]) / 2
    
    # Bin the events
    counts = fit(Histogram, event_times, bin_edges).weights
    
    # Calculate errors if requested
    count_errors = if show_errors
        if err_method == :poisson
            sqrt.(max.(counts, 1))  # Poisson errors, avoid sqrt(0)
        else  # gaussian
            sqrt.(counts)
        end
    else
        nothing
    end

    # Basic plot settings
    title --> "Event List Light Curve"
    xlabel --> "Time"
    ylabel --> "Counts per bin"
    grid --> true
    minorgrid --> true
    legend --> :topright

    # Axis limits handling
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits
            if !isnothing(xmin) || !isnothing(xmax)
                xlims --> (
                    isnothing(xmin) ? plot_tstart : xmin,
                    isnothing(xmax) ? plot_tstop : xmax,
                )
            end
            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? 0 : ymin,
                    isnothing(ymax) ? maximum(counts) * 1.1 : ymax,
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? plot_tstart : xmin,
                isnothing(xmax) ? plot_tstop : xmax,
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    end

    # Handle GTI/BTI visualization
    if show_btis || show_gtis
        effective_gtis = nothing

        # Priority: explicit gtis > gti_file > eventlist.meta.gti
        if !isnothing(gtis)
            effective_gtis = gtis
        elseif !isnothing(gti_file)
            try
                effective_gtis, _ = read_gti_from_fits(gti_file, gti_hdu_candidates=[gti_hdu])
            catch e
                @warn "Could not load GTIs from file: $e"
            end
        elseif has_gti(el)
            effective_gtis = gti(el)
        end

        if !isnothing(effective_gtis)
            y_min = 0
            y_max = maximum(counts) * 1.1

            if show_gtis
                # Plot GTI intervals as green rectangles
                for i in 1:size(effective_gtis, 1)
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]
                    
                    # Only show GTI intervals that overlap with plot range
                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start_clipped = max(gti_start, plot_tstart)
                        gti_stop_clipped = min(gti_stop, plot_tstop)
                        
                        @series begin
                            seriestype := :shape
                            fillcolor := :green
                            fillalpha := gti_alpha
                            linecolor := :green
                            linewidth := 0.5
                            label := i == 1 ? "Good Time Intervals" : ""
                            
                            # Rectangle coordinates
                            x_coords = [gti_start_clipped, gti_stop_clipped, gti_stop_clipped, gti_start_clipped, gti_start_clipped]
                            y_coords = [y_min, y_min, y_max, y_max, y_min]
                            (x_coords, y_coords)
                        end
                    end
                end
            end

            if show_btis
                # Calculate BTI intervals (gaps between GTIs)
                btis = []
                
                # Add BTI before first GTI if needed
                if effective_gtis[1, 1] > plot_tstart
                    push!(btis, [plot_tstart, effective_gtis[1, 1]])
                end
                
                # Add BTIs between GTI intervals
                for i in 1:(size(effective_gtis, 1) - 1)
                    gap_start = effective_gtis[i, 2]
                    gap_end = effective_gtis[i+1, 1]
                    if gap_end > gap_start && gap_end - gap_start > gap_threshold
                        push!(btis, [gap_start, gap_end])
                    end
                end
                
                # Add BTI after last GTI if needed
                if effective_gtis[end, 2] < plot_tstop
                    push!(btis, [effective_gtis[end, 2], plot_tstop])
                end

                # Plot BTI intervals as red rectangles
                for (i, bti) in enumerate(btis)
                    bti_start, bti_stop = bti[1], bti[2]
                    
                    if bti_stop > bti_start
                        @series begin
                            seriestype := :shape
                            fillcolor := :red
                            fillalpha := bti_alpha
                            linecolor := :red
                            linewidth := 0.5
                            label := i == 1 ? "Bad Time Intervals" : ""
                            
                            # Rectangle coordinates
                            x_coords = [bti_start, bti_stop, bti_stop, bti_start, bti_start]
                            y_coords = [y_min, y_min, y_max, y_max, y_min]
                            (x_coords, y_coords)
                        end
                    end
                end
            end
        else
            @debug "No GTI data found for visualization"
        end
    end

    # Main histogram plot
    if show_errors && !isnothing(count_errors)
        @series begin
            seriestype := :scatter
            marker := :circle
            markersize := 3
            markercolor := :blue
            markerstrokewidth := 0.5
            yerror := count_errors
            errorbar_color := :black
            color := :blue
            label := "Event Histogram with $(err_method) errors ($(bin_size)s bins)"
            (bin_centers, counts)
        end
    else
        @series begin
            seriestype := :steppost
            linewidth := 2
            color := :blue
            label := "Event Histogram ($(bin_size)s bins, $(length(event_times)) events)"
            
            # Create step plot coordinates
            step_x = Float64[]
            step_y = Float64[]
            
            for i in eachindex(bin_centers)
                bin_left = bin_edges[i]
                bin_right = bin_edges[i+1]
                
                push!(step_x, bin_left)
                push!(step_y, counts[i])
                push!(step_x, bin_right)
                push!(step_y, counts[i])
            end
            
            (step_x, step_y)
        end
    end
    
    # Add summary info as annotation
    @series begin
        seriestype := :scatter
        marker := :none
        label := ""
        annotations := [(plot_tstart + 0.05 * (plot_tstop - plot_tstart), 
                        maximum(counts) * 0.9,
                        text("$(length(event_times)) events\nBin size: $(bin_size)s", 
                             :left, 8, :gray))]
        ([plot_tstart], [maximum(counts) * 0.9])
    end
end

"""
    plot(lc::LightCurve{T}; kwargs...)

Plot a pre-computed light curve with optional properties, GTIs, and BTIs.

# Arguments
- `lc::LightCurve{T}`: Light curve object containing binned time series data

# Keywords
- `show_errors=false`: Display error bars if available
- `show_properties=false`: Show additional properties on secondary y-axis
- `property_name=:mean_energy`: Which property to display (if `show_properties=true`)
- `show_btis=false`: Show Bad Time Intervals as red shaded regions
- `show_bti=false`: Alias for `show_btis`
- `show_gtis=false`: Show Good Time Intervals as green shaded regions
- `show_gti=false`: Alias for `show_gtis`
- `gtis=nothing`: GTI matrix to use (overrides metadata GTIs)
- `gti_file=nothing`: FITS file containing GTI extension
- `gti_hdu="GTI"`: HDU name for GTI data
- `bti_alpha=0.3`: Transparency for BTI shading
- `gti_alpha=0.2`: Transparency for GTI shading
- `axis_limits=nothing`: Plot limits as `[xmin, xmax]` or `[xmin, xmax, ymin, ymax]`

# Returns
- `Tuple{Vector, Vector}`: Time bins and corresponding count rates

# Examples
```julia
# Basic light curve plot
plot(lightcurve)

# With error bars and mean energy overlay
plot(lightcurve, show_errors=true, show_properties=true, property_name=:mean_energy)

# Show GTIs with custom transparency
plot(lightcurve, show_gtis=true, gti_alpha=0.4)
```

# Notes
- Errors are calculated on-demand if not already present
- Properties must exist in `lc.properties` to be displayed
- GTI priority: explicit `gtis` > `gti_file` > `lc.metadata.extra["gti_bounds"]`
"""
@recipe function f(
    lc::LightCurve{T};
    tstart = nothing,
    tstop = nothing,
    show_errors = false,
    show_properties = false,
    property_name = :mean_energy,
    show_btis = false,
    show_bti = false,
    show_gtis = false,
    show_gti = false,
    gtis = nothing,
    gti_file = nothing,
    gti_hdu = "GTI",
    bti_alpha = 0.3,
    gti_alpha = 0.2,
    axis_limits = nothing,
) where {T}

    show_gtis = show_gtis || show_gti
    show_btis = show_btis || show_bti

    if show_errors && isnothing(lc.count_error)
        calculate_errors!(lc)
    end

    # Apply time filtering if specified
    time_mask = trues(length(lc.time))
    
    if !isnothing(tstart)
        time_mask .&= (lc.time .>= tstart)
    end
    
    if !isnothing(tstop)
        time_mask .&= (lc.time .<= tstop)
    end
    
    # Check if any data remains after filtering
    if !any(time_mask)
        @warn "No data points in specified time range [$tstart, $tstop]"
        # Create empty plot
        title --> "Light Curve (No Data in Range)"
        xlabel --> "Time (s)"
        ylabel --> "Counts"
        @series begin
            seriestype := :scatter
            marker := :none
            label := ""
            ([0], [0])
        end
        return
    end
    
    # Filter the data
    filtered_time = lc.time[time_mask]
    filtered_counts = lc.counts[time_mask]
    filtered_errors = isnothing(lc.count_error) ? nothing : lc.count_error[time_mask]
    
    # Convert to matrix format
    lc_matrix = if isnothing(filtered_errors)
        hcat(filtered_time, filtered_counts, zeros(length(filtered_time)))
    else
        hcat(filtered_time, filtered_counts, filtered_errors)
    end

    title --> "Light Curve"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    minorgrid --> true
    legend --> :bottomright

    # Axis limits handling
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits

            if !isnothing(xmin) || !isnothing(xmax)
                xlims --> (
                    isnothing(xmin) ? minimum(lc_matrix[:, 1]) : xmin,
                    isnothing(xmax) ? maximum(lc_matrix[:, 1]) : xmax,
                )
            end

            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(lc_matrix[:, 2]) : ymin,
                    isnothing(ymax) ? maximum(lc_matrix[:, 2]) : ymax,
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? minimum(lc_matrix[:, 1]) : xmin,
                isnothing(xmax) ? maximum(lc_matrix[:, 1]) : xmax,
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    else
        # Set automatic limits based on filtered data
        xlims --> (minimum(filtered_time), maximum(filtered_time))
    end

    # Determine plot time range (use filtered data or tstart/tstop)
    plot_tstart = isnothing(tstart) ? minimum(filtered_time) : tstart
    plot_tstop = isnothing(tstop) ? maximum(filtered_time) : tstop

    # Handle GTI/BTI visualization
    if show_btis || show_gtis
        effective_gtis = nothing

        if !isnothing(gtis)
            effective_gtis = gtis
        elseif !isnothing(gti_file)
            try
                effective_gtis = load_gtis(gti_file, gti_hdu)
            catch e
                @warn "Could not load GTIs from file: $e"
            end
        elseif haskey(lc.metadata.extra, "gti") && !isnothing(lc.metadata.extra["gti"])
            effective_gtis = lc.metadata.extra["gti"]
        end

        if !isnothing(effective_gtis)
            y_min, y_max = extrema(lc_matrix[:, 2])

            if show_gtis
                n_gtis = size(effective_gtis, 1)
                gti_x = Vector{Float64}(undef, n_gtis * 6)
                gti_y = Vector{Float64}(undef, n_gtis * 6)
                gti_idx = 0

                @inbounds for i = 1:n_gtis
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]

                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start_clipped = max(gti_start, plot_tstart)
                        gti_stop_clipped = min(gti_stop, plot_tstop)

                        base_idx = gti_idx * 6
                        gti_x[base_idx+1] = gti_start_clipped
                        gti_x[base_idx+2] = gti_stop_clipped
                        gti_x[base_idx+3] = gti_stop_clipped
                        gti_x[base_idx+4] = gti_start_clipped
                        gti_x[base_idx+5] = gti_start_clipped
                        gti_x[base_idx+6] = NaN

                        gti_y[base_idx+1] = y_min
                        gti_y[base_idx+2] = y_min
                        gti_y[base_idx+3] = y_max
                        gti_y[base_idx+4] = y_max
                        gti_y[base_idx+5] = y_min
                        gti_y[base_idx+6] = NaN

                        gti_idx += 1
                    end
                end

                if gti_idx > 0
                    resize!(gti_x, gti_idx * 6)
                    resize!(gti_y, gti_idx * 6)

                    @series begin
                        seriestype := :shape
                        fillcolor := :green
                        fillalpha := gti_alpha
                        linecolor := :green
                        linewidth := 0.5
                        label := "Good Time Intervals"
                        gti_x, gti_y
                    end
                end
            end

            if show_btis
                btis = get_btis(effective_gtis, plot_tstart, plot_tstop)

                if !isempty(btis)
                    n_btis = size(btis, 1)
                    bti_x = Vector{Float64}(undef, n_btis * 6)
                    bti_y = Vector{Float64}(undef, n_btis * 6)

                    @inbounds for i = 1:n_btis
                        bti_start, bti_stop = btis[i, 1], btis[i, 2]

                        base_idx = (i - 1) * 6
                        bti_x[base_idx+1] = bti_start
                        bti_x[base_idx+2] = bti_stop
                        bti_x[base_idx+3] = bti_stop
                        bti_x[base_idx+4] = bti_start
                        bti_x[base_idx+5] = bti_start
                        bti_x[base_idx+6] = NaN

                        bti_y[base_idx+1] = y_min
                        bti_y[base_idx+2] = y_min
                        bti_y[base_idx+3] = y_max
                        bti_y[base_idx+4] = y_max
                        bti_y[base_idx+5] = y_min
                        bti_y[base_idx+6] = NaN
                    end

                    @series begin
                        seriestype := :shape
                        fillcolor := :red
                        fillalpha := bti_alpha
                        linecolor := :red
                        linewidth := 0.5
                        label := "Bad Time Intervals"
                        bti_x, bti_y
                    end
                end
            end
        end
    end

    # Handle properties display (filter properties too)
    if show_properties
        prop_idx = findfirst(p -> p.name == property_name, lc.properties)
        if !isnothing(prop_idx)
            prop = lc.properties[prop_idx]
            filtered_prop_values = prop.values[time_mask]
            
            prop_matrix = hcat(filtered_time, filtered_prop_values)

            @series begin
                yaxis := :right
                ylabel := "$(prop.name) ($(prop.unit))"
                seriestype := :line
                color := :red
                linewidth := 1.5
                label := String(prop.name)
                prop_matrix[:, 1], prop_matrix[:, 2]
            end
        end
    end

    # Main light curve plot
    if show_errors && !isnothing(filtered_errors)
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:, 3]
        errorbar_color --> :black
        color --> :blue
        label --> "Light Curve with $(lc.err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Light Curve (bin size: $(lc.metadata.bin_size)s)"
    end

    return lc_matrix[:, 1], lc_matrix[:, 2]
end

"""
    plot(lc_segments::Vector{<:LightCurve}; kwargs...)

Plot multiple light curve segments with optional segment boundaries and individual coloring.

# Arguments
- `lc_segments::Vector{<:LightCurve}`: Vector of light curve segments to plot

# Keywords
- `show_errors=false`: Display error bars for each segment
- `show_segment_boundaries=true`: Show vertical dashed lines at segment boundaries
- `segment_colors=nothing`: Custom colors for each segment (defaults to standard palette)
- `axis_limits=nothing`: Plot limits as `[xmin, xmax]` or `[xmin, xmax, ymin, ymax]`

# Returns
- Multiple plot series, one per segment

# Examples
```julia
# Basic segmented light curve
plot(segments)

# With custom colors and no boundaries
plot(segments, segment_colors=[:red, :blue, :green], show_segment_boundaries=false)

# With error bars and custom limits
plot(segments, show_errors=true, axis_limits=[0, 1000, 0, 50])
```

# Notes
- Default color palette cycles through 8 colors for segments
- Segment boundaries are drawn at the start of each segment after the first
- Each segment can have different bin sizes and error methods
"""
@recipe function f(
    lc_segments::Vector{<:LightCurve};
    show_errors = false,
    show_segment_boundaries = true,
    segment_colors = nothing,
    axis_limits = nothing,
)

    title --> "Segmented Light Curve"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    minorgrid --> true
    legend --> :bottomright

    # Axis limits handling
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits

            # Get overall data bounds
            all_times = vcat([lc.time for lc in lc_segments]...)
            all_counts = vcat([lc.counts for lc in lc_segments]...)

            if !isnothing(xmin) || !isnothing(xmax)
                xlims --> (
                    isnothing(xmin) ? minimum(all_times) : xmin,
                    isnothing(xmax) ? maximum(all_times) : xmax,
                )
            end

            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(all_counts) : ymin,
                    isnothing(ymax) ? maximum(all_counts) : ymax,
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            all_times = vcat([lc.time for lc in lc_segments]...)
            xlims --> (
                isnothing(xmin) ? minimum(all_times) : xmin,
                isnothing(xmax) ? maximum(all_times) : xmax,
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    end

    default_colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
    colors = isnothing(segment_colors) ? default_colors : segment_colors
    n_colors = length(colors)

    boundaries = Vector{Float64}()

    for (i, lc) in enumerate(lc_segments)
        color = colors[((i-1)%n_colors)+1]

        if show_errors && isnothing(lc.count_error)
            calculate_errors!(lc)
        end

        lc_matrix = hcat(lc.time, lc.counts, lc.count_error)

        @series begin
            if show_errors
                seriestype := :scatter
                marker := :circle
                markersize := 3
                markerstrokewidth := 0.5
                yerror := lc_matrix[:, 3]
                errorbar_color := color
            else
                seriestype := :steppost
                linewidth := 1.5
            end
            color := color
            label := "Segment $i"

            lc_matrix[:, 1], lc_matrix[:, 2]
        end

        if show_segment_boundaries && i > 1
            push!(boundaries, minimum(lc.time))
        end
    end

    if show_segment_boundaries && !isempty(boundaries)
        @series begin
            seriestype := :vline
            color := :black
            linestyle := :dash
            linewidth := 1
            alpha := 0.7
            label := "Segment boundaries"
            boundaries
        end
    end
end
"""
    plot(lc::LightCurve{T}, new_binsize::Real; kwargs...)

Plot a light curve after rebinning it to a new bin size.

# Arguments
- `lc::LightCurve{T}`: Original light curve to rebin and plot
- `new_binsize::Real`: New bin size in seconds (must be larger than current bin size)

# Keywords
- `show_errors=false`: Display error bars using rebinned error estimates
- `show_properties=false`: Show additional properties on secondary y-axis
- `property_name=:mean_energy`: Which property to display (if `show_properties=true`)
- `show_btis=false`: Show Bad Time Intervals as red shaded regions
- `show_bti=false`: Alias for `show_btis`
- `show_gtis=false`: Show Good Time Intervals as green shaded regions
- `show_gti=false`: Alias for `show_gtis`
- `gtis=nothing`: GTI matrix to use (overrides metadata GTIs)
- `gti_file=nothing`: FITS file containing GTI extension
- `gti_hdu="GTI"`: HDU name for GTI data
- `bti_alpha=0.3`: Transparency for BTI shading
- `gti_alpha=0.2`: Transparency for GTI shading
- `axis_limits=nothing`: Plot limits as `[xmin, xmax]` or `[xmin, xmax, ymin, ymax]`
- `show_original=false`: Overlay original light curve for comparison
- `original_alpha=0.3`: Transparency for original light curve overlay

# Returns
- `Tuple{Vector, Vector}`: Rebinned time bins and corresponding count rates

# Examples
```julia
# Basic rebinned light curve
plot(lc, 100.0)  # Rebin to 100s

# With error bars and GTIs
plot(lc, 50.0, show_errors=true, show_gtis=true)

# Compare with original
plot(lc, 200.0, show_original=true, original_alpha=0.4)

# With properties overlay
plot(lc, 30.0, show_properties=true, property_name=:mean_energy)
```

# Notes
- The new bin size must be larger than the current bin size
- Original light curve data is preserved; only the plot shows rebinned data
- Error bars are recalculated for the rebinned light curve
- GTI/BTI regions are preserved from the original light curve metadata
"""
@recipe function f(
    lc::LightCurve{T},
    new_binsize::Real;
    show_errors = false,
    show_properties = false,
    property_name = :mean_energy,
    show_btis = false,
    show_bti = false,
    show_gtis = false,
    show_gti = false,
    gtis = nothing,
    gti_file = nothing,
    gti_hdu = "GTI",
    bti_alpha = 0.3,
    gti_alpha = 0.2,
    axis_limits = nothing,
    show_original = false,
    original_alpha = 0.3,
) where {T}

    # Validate new bin size
    new_binsize <= lc.metadata.bin_size && throw(
        ArgumentError(
            "New bin size ($new_binsize s) must be larger than current bin size ($(lc.metadata.bin_size) s)",
        ),
    )

    # Create rebinned light curve
    rebinned_lc = rebin(lc, new_binsize)

    show_gtis = show_gtis || show_gti
    show_btis = show_btis || show_bti

    if show_errors && isnothing(rebinned_lc.count_error)
        calculate_errors!(rebinned_lc)
    end

    # Convert to matrix format
    lc_matrix = hcat(rebinned_lc.time, rebinned_lc.counts, rebinned_lc.count_error)

    title --> "Rebinned Light Curve ($(lc.metadata.bin_size)s → $(new_binsize)s)"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    minorgrid --> true
    legend --> :bottomright

    # Axis limits handling
    if !isnothing(axis_limits)
        if length(axis_limits) == 4
            xmin, xmax, ymin, ymax = axis_limits

            if !isnothing(xmin) || !isnothing(xmax)
                xlims --> (
                    isnothing(xmin) ? minimum(lc_matrix[:, 1]) : xmin,
                    isnothing(xmax) ? maximum(lc_matrix[:, 1]) : xmax,
                )
            end

            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(lc_matrix[:, 2]) : ymin,
                    isnothing(ymax) ? maximum(lc_matrix[:, 2]) : ymax,
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? minimum(lc_matrix[:, 1]) : xmin,
                isnothing(xmax) ? maximum(lc_matrix[:, 1]) : xmax,
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    end

    plot_tstart, plot_tstop = rebinned_lc.metadata.time_range

    # Show original light curve for comparison
    if show_original
        original_matrix = hcat(lc.time, lc.counts, lc.count_error)

        @series begin
            seriestype := :steppost
            linewidth := 1
            color := :gray
            alpha := original_alpha
            label := "Original ($(lc.metadata.bin_size)s bins)"
            original_matrix[:, 1], original_matrix[:, 2]
        end
    end

    # Handle GTI/BTI visualization (use original metadata)
    if show_btis || show_gtis
        effective_gtis = nothing

        if !isnothing(gtis)
            effective_gtis = gtis
        elseif !isnothing(gti_file)
            try
                effective_gtis = load_gtis(gti_file, gti_hdu)
            catch e
                @warn "Could not load GTIs from file: $e"
            end
        elseif haskey(lc.metadata.extra, "gti_applied") &&
               haskey(lc.metadata.extra, "gti_bounds")
            gti_bounds = lc.metadata.extra["gti_bounds"]
            effective_gtis = reshape(gti_bounds, 1, 2)
        end

        if !isnothing(effective_gtis)
            y_min, y_max = extrema(lc_matrix[:, 2])

            if show_gtis
                n_gtis = size(effective_gtis, 1)
                gti_x = Vector{Float64}(undef, n_gtis * 6)
                gti_y = Vector{Float64}(undef, n_gtis * 6)
                gti_idx = 0

                @inbounds for i = 1:n_gtis
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]

                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start = max(gti_start, plot_tstart)
                        gti_stop = min(gti_stop, plot_tstop)

                        base_idx = gti_idx * 6
                        gti_x[base_idx+1] = gti_start
                        gti_x[base_idx+2] = gti_stop
                        gti_x[base_idx+3] = gti_stop
                        gti_x[base_idx+4] = gti_start
                        gti_x[base_idx+5] = gti_start
                        gti_x[base_idx+6] = NaN

                        gti_y[base_idx+1] = y_min
                        gti_y[base_idx+2] = y_min
                        gti_y[base_idx+3] = y_max
                        gti_y[base_idx+4] = y_max
                        gti_y[base_idx+5] = y_min
                        gti_y[base_idx+6] = NaN

                        gti_idx += 1
                    end
                end

                if gti_idx > 0
                    resize!(gti_x, gti_idx * 6)
                    resize!(gti_y, gti_idx * 6)

                    @series begin
                        seriestype := :shape
                        fillcolor := :green
                        fillalpha := gti_alpha
                        linecolor := :green
                        linewidth := 0.5
                        label := "Good Time Intervals"
                        gti_x, gti_y
                    end
                end
            end

            if show_btis
                btis = get_btis(effective_gtis, plot_tstart, plot_tstop)

                if !isempty(btis)
                    n_btis = size(btis, 1)
                    bti_x = Vector{Float64}(undef, n_btis * 6)
                    bti_y = Vector{Float64}(undef, n_btis * 6)

                    @inbounds for i = 1:n_btis
                        bti_start, bti_stop = btis[i, 1], btis[i, 2]

                        base_idx = (i - 1) * 6
                        bti_x[base_idx+1] = bti_start
                        bti_x[base_idx+2] = bti_stop
                        bti_x[base_idx+3] = bti_stop
                        bti_x[base_idx+4] = bti_start
                        bti_x[base_idx+5] = bti_start
                        bti_x[base_idx+6] = NaN

                        bti_y[base_idx+1] = y_min
                        bti_y[base_idx+2] = y_min
                        bti_y[base_idx+3] = y_max
                        bti_y[base_idx+4] = y_max
                        bti_y[base_idx+5] = y_min
                        bti_y[base_idx+6] = NaN
                    end

                    @series begin
                        seriestype := :shape
                        fillcolor := :red
                        fillalpha := bti_alpha
                        linecolor := :red
                        linewidth := 0.5
                        label := "Bad Time Intervals"
                        bti_x, bti_y
                    end
                end
            end
        end
    end

    # Handle properties display
    if show_properties
        prop_idx = findfirst(p -> p.name == property_name, rebinned_lc.properties)
        if !isnothing(prop_idx)
            prop = rebinned_lc.properties[prop_idx]

            prop_matrix = hcat(rebinned_lc.time, prop.values)

            @series begin
                yaxis := :right
                ylabel := "$(prop.name) ($(prop.unit))"
                seriestype := :line
                color := :red
                linewidth := 1.5
                label := String(prop.name)
                prop_matrix[:, 1], prop_matrix[:, 2]
            end
        end
    end

    if show_errors
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:, 3]
        errorbar_color --> :black
        color --> :blue
        label --> "Rebinned LC ($(new_binsize)s) with $(rebinned_lc.err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Rebinned LC ($(new_binsize)s bins)"
    end

    return lc_matrix[:, 1], lc_matrix[:, 2]
end
