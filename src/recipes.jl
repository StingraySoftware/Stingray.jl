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
@recipe function f(el::EventList{T}, bin_size::Real=1.0;
                  tstart=nothing,
                  tstop=nothing,
                  energy_filter=nothing,
                  show_errors=false,
                  show_btis=false,
                  show_bti=false,
                  show_gtis=false,
                  show_gti=false,
                  gtis=nothing,
                  gti_file=nothing,
                  gti_hdu="GTI",
                  bti_alpha=0.3,
                  gti_alpha=0.2,
                  gap_threshold=10.0,
                  axis_limits=nothing,
                  err_method=:poisson) where T

    isempty(el.times) && error("EventList is empty")
    
    show_gtis = show_gtis || show_gti
    show_btis = show_btis || show_bti
    
    # Create light curve
    lc = create_lightcurve(el, bin_size;
                          tstart=tstart,
                          tstop=tstop,
                          energy_filter=energy_filter,
                          err_method=err_method)

    calculate_errors!(lc)

    # Convert to matrix format
    lc_matrix = hcat(lc.time, lc.counts, lc.count_error)
    
    # Basic plot settings
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
                    isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                    isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
                )
            end
            
            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(lc_matrix[:,2]) : ymin,
                    isnothing(ymax) ? maximum(lc_matrix[:,2]) : ymax
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    end

    # Determine time range
    plot_tstart = isnothing(tstart) ? lc.time[1] : tstart
    plot_tstop = isnothing(tstop) ? lc.time[end] : tstop
    
    # Handle GTI/BTI visualization
    if show_btis || show_gtis
        effective_gtis = nothing
        
        # Priority: explicit gtis > gti_file > eventlist.meta.gti
        if !isnothing(gtis)
            effective_gtis = gtis
        elseif !isnothing(gti_file)
            try
                effective_gtis = load_gtis(gti_file, gti_hdu)
            catch e
                @warn "Could not load GTIs from file: $e"
            end
        elseif !isnothing(el.meta.gti)
            effective_gtis = el.meta.gti
        end
        
        if !isnothing(effective_gtis)
            y_min, y_max = extrema(lc_matrix[:,2])
            
            # Pre-allocate arrays for shapes
            if show_gtis
                n_gtis = size(effective_gtis, 1)
                gti_x = Vector{Float64}(undef, n_gtis * 6)
                gti_y = Vector{Float64}(undef, n_gtis * 6)
                gti_idx = 0
                
                @inbounds for i in 1:n_gtis
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]
                    
                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start = max(gti_start, plot_tstart)
                        gti_stop = min(gti_stop, plot_tstop)
                        
                        # Rectangle vertices
                        base_idx = gti_idx * 6
                        gti_x[base_idx + 1] = gti_start
                        gti_x[base_idx + 2] = gti_stop
                        gti_x[base_idx + 3] = gti_stop
                        gti_x[base_idx + 4] = gti_start
                        gti_x[base_idx + 5] = gti_start
                        gti_x[base_idx + 6] = NaN
                        
                        gti_y[base_idx + 1] = y_min
                        gti_y[base_idx + 2] = y_min
                        gti_y[base_idx + 3] = y_max
                        gti_y[base_idx + 4] = y_max
                        gti_y[base_idx + 5] = y_min
                        gti_y[base_idx + 6] = NaN
                        
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
                    
                    @inbounds for i in 1:n_btis
                        bti_start, bti_stop = btis[i, 1], btis[i, 2]
                        
                        base_idx = (i - 1) * 6
                        bti_x[base_idx + 1] = bti_start
                        bti_x[base_idx + 2] = bti_stop
                        bti_x[base_idx + 3] = bti_stop
                        bti_x[base_idx + 4] = bti_start
                        bti_x[base_idx + 5] = bti_start
                        bti_x[base_idx + 6] = NaN
                        
                        bti_y[base_idx + 1] = y_min
                        bti_y[base_idx + 2] = y_min
                        bti_y[base_idx + 3] = y_max
                        bti_y[base_idx + 4] = y_max
                        bti_y[base_idx + 5] = y_min
                        bti_y[base_idx + 6] = NaN
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

    # Main light curve plot
    if show_errors
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:,3]
        errorbar_color --> :black
        color --> :blue
        label --> "Light Curve with $(err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Light Curve (bin size: $bin_size s)"
    end
    
    return lc_matrix[:,1], lc_matrix[:,2]
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
@recipe function f(lc::LightCurve{T};
                  show_errors=false,
                  show_properties=false,
                  property_name=:mean_energy,
                  show_btis=false,
                  show_bti=false,
                  show_gtis=false,
                  show_gti=false,
                  gtis=nothing,
                  gti_file=nothing,
                  gti_hdu="GTI",
                  bti_alpha=0.3,
                  gti_alpha=0.2,
                  axis_limits=nothing) where T

    show_gtis = show_gtis || show_gti
    show_btis = show_btis || show_bti

    if show_errors && isnothing(lc.count_error)
        calculate_errors!(lc)
    end

    # Convert to matrix format
    lc_matrix = hcat(lc.time, lc.counts, lc.count_error)
    
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
                    isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                    isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
                )
            end
            
            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(lc_matrix[:,2]) : ymin,
                    isnothing(ymax) ? maximum(lc_matrix[:,2]) : ymax
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
            )
        else
            @warn "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]"
        end
    end

    plot_tstart, plot_tstop = lc.metadata.time_range
    
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
        elseif haskey(lc.metadata.extra, "gti_applied") && haskey(lc.metadata.extra, "gti_bounds")
            gti_bounds = lc.metadata.extra["gti_bounds"]
            effective_gtis = reshape(gti_bounds, 1, 2)
        end
        
        if !isnothing(effective_gtis)
            y_min, y_max = extrema(lc_matrix[:,2])
            
            if show_gtis
                n_gtis = size(effective_gtis, 1)
                gti_x = Vector{Float64}(undef, n_gtis * 6)
                gti_y = Vector{Float64}(undef, n_gtis * 6)
                gti_idx = 0
                
                @inbounds for i in 1:n_gtis
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]
                    
                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start = max(gti_start, plot_tstart)
                        gti_stop = min(gti_stop, plot_tstop)
                        
                        base_idx = gti_idx * 6
                        gti_x[base_idx + 1] = gti_start
                        gti_x[base_idx + 2] = gti_stop
                        gti_x[base_idx + 3] = gti_stop
                        gti_x[base_idx + 4] = gti_start
                        gti_x[base_idx + 5] = gti_start
                        gti_x[base_idx + 6] = NaN
                        
                        gti_y[base_idx + 1] = y_min
                        gti_y[base_idx + 2] = y_min
                        gti_y[base_idx + 3] = y_max
                        gti_y[base_idx + 4] = y_max
                        gti_y[base_idx + 5] = y_min
                        gti_y[base_idx + 6] = NaN
                        
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
                    
                    @inbounds for i in 1:n_btis
                        bti_start, bti_stop = btis[i, 1], btis[i, 2]
                        
                        base_idx = (i - 1) * 6
                        bti_x[base_idx + 1] = bti_start
                        bti_x[base_idx + 2] = bti_stop
                        bti_x[base_idx + 3] = bti_stop
                        bti_x[base_idx + 4] = bti_start
                        bti_x[base_idx + 5] = bti_start
                        bti_x[base_idx + 6] = NaN
                        
                        bti_y[base_idx + 1] = y_min
                        bti_y[base_idx + 2] = y_min
                        bti_y[base_idx + 3] = y_max
                        bti_y[base_idx + 4] = y_max
                        bti_y[base_idx + 5] = y_min
                        bti_y[base_idx + 6] = NaN
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
        prop_idx = findfirst(p -> p.name == property_name, lc.properties)
        if !isnothing(prop_idx)
            prop = lc.properties[prop_idx]
            
            prop_matrix = hcat(lc.time, prop.values)
            
            @series begin
                yaxis := :right
                ylabel := "$(prop.name) ($(prop.unit))"
                seriestype := :line
                color := :red
                linewidth := 1.5
                label := String(prop.name)
                prop_matrix[:,1], prop_matrix[:,2]
            end
        end
    end

    # Main light curve plot
    if show_errors
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:,3]
        errorbar_color --> :black
        color --> :blue
        label --> "Light Curve with $(lc.err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Light Curve (bin size: $(lc.metadata.bin_size)s)"
    end

    return lc_matrix[:,1], lc_matrix[:,2]
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
@recipe function f(lc_segments::Vector{<:LightCurve};
                  show_errors=false,
                  show_segment_boundaries=true,
                  segment_colors=nothing,
                  axis_limits=nothing)

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
                    isnothing(xmax) ? maximum(all_times) : xmax
                )
            end
            
            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(all_counts) : ymin,
                    isnothing(ymax) ? maximum(all_counts) : ymax
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            all_times = vcat([lc.time for lc in lc_segments]...)
            xlims --> (
                isnothing(xmin) ? minimum(all_times) : xmin,
                isnothing(xmax) ? maximum(all_times) : xmax
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
        color = colors[((i-1) % n_colors) + 1]
        
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
                yerror := lc_matrix[:,3]
                errorbar_color := color
            else
                seriestype := :steppost
                linewidth := 1.5
            end
            color := color
            label := "Segment $i"
            
            lc_matrix[:,1], lc_matrix[:,2]
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
@recipe function f(lc::LightCurve{T}, new_binsize::Real;
                  show_errors=false,
                  show_properties=false,
                  property_name=:mean_energy,
                  show_btis=false,
                  show_bti=false,
                  show_gtis=false,
                  show_gti=false,
                  gtis=nothing,
                  gti_file=nothing,
                  gti_hdu="GTI",
                  bti_alpha=0.3,
                  gti_alpha=0.2,
                  axis_limits=nothing,
                  show_original=false,
                  original_alpha=0.3) where T

    # Validate new bin size
    new_binsize <= lc.metadata.bin_size && throw(ArgumentError(
        "New bin size ($new_binsize s) must be larger than current bin size ($(lc.metadata.bin_size) s)"
    ))
    
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
                    isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                    isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
                )
            end
            
            if !isnothing(ymin) || !isnothing(ymax)
                ylims --> (
                    isnothing(ymin) ? minimum(lc_matrix[:,2]) : ymin,
                    isnothing(ymax) ? maximum(lc_matrix[:,2]) : ymax
                )
            end
        elseif length(axis_limits) == 2
            xmin, xmax = axis_limits
            xlims --> (
                isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin,
                isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax
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
            original_matrix[:,1], original_matrix[:,2]
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
        elseif haskey(lc.metadata.extra, "gti_applied") && haskey(lc.metadata.extra, "gti_bounds")
            gti_bounds = lc.metadata.extra["gti_bounds"]
            effective_gtis = reshape(gti_bounds, 1, 2)
        end
        
        if !isnothing(effective_gtis)
            y_min, y_max = extrema(lc_matrix[:,2])
            
            if show_gtis
                n_gtis = size(effective_gtis, 1)
                gti_x = Vector{Float64}(undef, n_gtis * 6)
                gti_y = Vector{Float64}(undef, n_gtis * 6)
                gti_idx = 0
                
                @inbounds for i in 1:n_gtis
                    gti_start, gti_stop = effective_gtis[i, 1], effective_gtis[i, 2]
                    
                    if gti_stop >= plot_tstart && gti_start <= plot_tstop
                        gti_start = max(gti_start, plot_tstart)
                        gti_stop = min(gti_stop, plot_tstop)
                        
                        base_idx = gti_idx * 6
                        gti_x[base_idx + 1] = gti_start
                        gti_x[base_idx + 2] = gti_stop
                        gti_x[base_idx + 3] = gti_stop
                        gti_x[base_idx + 4] = gti_start
                        gti_x[base_idx + 5] = gti_start
                        gti_x[base_idx + 6] = NaN
                        
                        gti_y[base_idx + 1] = y_min
                        gti_y[base_idx + 2] = y_min
                        gti_y[base_idx + 3] = y_max
                        gti_y[base_idx + 4] = y_max
                        gti_y[base_idx + 5] = y_min
                        gti_y[base_idx + 6] = NaN
                        
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
                    
                    @inbounds for i in 1:n_btis
                        bti_start, bti_stop = btis[i, 1], btis[i, 2]
                        
                        base_idx = (i - 1) * 6
                        bti_x[base_idx + 1] = bti_start
                        bti_x[base_idx + 2] = bti_stop
                        bti_x[base_idx + 3] = bti_stop
                        bti_x[base_idx + 4] = bti_start
                        bti_x[base_idx + 5] = bti_start
                        bti_x[base_idx + 6] = NaN
                        
                        bti_y[base_idx + 1] = y_min
                        bti_y[base_idx + 2] = y_min
                        bti_y[base_idx + 3] = y_max
                        bti_y[base_idx + 4] = y_max
                        bti_y[base_idx + 5] = y_min
                        bti_y[base_idx + 6] = NaN
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
                prop_matrix[:,1], prop_matrix[:,2]
            end
        end
    end

    if show_errors
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:,3]
        errorbar_color --> :black
        color --> :blue
        label --> "Rebinned LC ($(new_binsize)s) with $(rebinned_lc.err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Rebinned LC ($(new_binsize)s bins)"
    end

    return lc_matrix[:,1], lc_matrix[:,2]
end