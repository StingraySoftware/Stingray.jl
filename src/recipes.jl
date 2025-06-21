@recipe function f(el::EventList{T}, bin_size::Real=1.0;
                  tstart=nothing,
                  tstop=nothing,
                  energy_filter=nothing,
                  show_errors=false,
                  show_gaps=false,
                  gap_threshold=10.0,
                  axis_limits=nothing,
                  err_method=:poisson) where T
    isempty(el.times) && error("EventList is empty")
    lc = create_lightcurve(el, bin_size;
                          tstart=tstart,
                          tstop=tstop,
                          energy_filter=energy_filter,
                          err_method=err_method)

    lc_matrix = hcat(lc.timebins, lc.counts, lc.count_error)

    title --> "Light Curve"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    minorgrid --> true
    legend --> :bottomright
    if !isnothing(axis_limits) && length(axis_limits) == 4
        xmin, xmax, ymin, ymax = axis_limits
        xlims --> (isnothing(xmin) ? minimum(lc_matrix[:,1]) : xmin, 
                  isnothing(xmax) ? maximum(lc_matrix[:,1]) : xmax)
        ylims --> (isnothing(ymin) ? minimum(lc_matrix[:,2]) : ymin, 
                  isnothing(ymax) ? maximum(lc_matrix[:,2]) : ymax)
    end

    if show_gaps
        if !isnothing(energy_filter) && !isnothing(el.energies)
            emin, emax = energy_filter
            energy_mask = (el.energies .>= emin) .& (el.energies .< emax)
            filtered_times = el.times[energy_mask]
        else
            filtered_times = el.times
        end
        
        sorted_times = sort(filtered_times)
        time_diffs = diff(sorted_times)
        median_diff = median(time_diffs)
        gap_indices = findall(x -> x > (gap_threshold * median_diff), time_diffs)

        if !isempty(gap_indices)
            @series begin
                seriestype := :vline
                color := :red
                linewidth := 1
                label := "Gaps"
                sorted_times[gap_indices]
            end
        end
    end
    if show_errors
        seriestype --> :scatter
        marker --> :circle
        markersize --> 3
        markercolor --> :blue
        markerstrokewidth --> 0.5
        yerror --> lc_matrix[:,3]  # Error column
        errorbar_color --> :black
        color --> :blue
        label --> "Light Curve with $(err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Light Curve (bin size: $bin_size)"
    end

    return lc_matrix[:,1], lc_matrix[:,2]  # time bins and counts
end

@recipe function f(lc::LightCurve{T};
                  show_errors=false,
                  show_properties=false,
                  property_name=:mean_energy) where T

    # Convert to matrix format
    lc_matrix = hcat(lc.timebins, lc.counts, lc.count_error)

    title --> "Light Curve"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    grid --> true
    minorgrid --> true
    legend --> :bottomright

    if show_properties
        prop_idx = findfirst(p -> p.name == property_name, lc.properties)
        if !isnothing(prop_idx)
            prop = lc.properties[prop_idx]
            
            # Property data as matrix
            prop_matrix = hcat(lc.timebins, prop.values)
            
            @series begin
                yaxis := :right
                ylabel := "$(prop.name) ($(prop.unit))"
                seriestype := :line
                color := :red
                linewidth := 1.5
                label := prop.name
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
        label --> "Light Curve with $(lc.err_method) errors"
    else
        seriestype --> :steppost
        linewidth --> 1.5
        color --> :blue
        label --> "Light Curve (bin size: $(lc.metadata.bin_size))"
    end

    return lc_matrix[:,1], lc_matrix[:,2]
end