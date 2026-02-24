using RecipesBase
using Statistics

@recipe function f(events::EventList, binsize::Real=1.0; 
                  energy_filter=nothing, 
                  show_gtis=false)
    
    # Defaults
    title  --> "EventList Lightcurve (bin: $(binsize)s)"
    xlabel --> "Time (s)"
    ylabel --> "Counts/bin"
    seriestype --> :steppre

    # 1. Filter times based on energy if requested
    times = events.times
    if !isnothing(energy_filter) && !isnothing(events.energies)
        mask = (events.energies .>= energy_filter[1]) .& (events.energies .< energy_filter[2])
        times = times[mask]
    end

    # 2. GTI Layer (Our "Art")
    if show_gtis && has_gti(events)
        # We loop through each interval in the GTI matrix
        gtis = events.meta.gti
        for i in 1:size(gtis, 1)
            @series begin
                seriestype := :vspan
                fillcolor := :green
                fillalpha := 0.15
                label := (i == 1 ? "GTI" : "") # Only label the first one for the legend
                primary := false
                [gtis[i, 1], gtis[i, 2]]
            end
        end
    end
    # 3. Innovation: Manual Binning
    # We create the bins and count how many events fall into each
    if isempty(times)
        return Float64[], Float64[]
    end
    
    t_start, t_end = minimum(times), maximum(times)
    edges = collect(t_start:binsize:(t_end + binsize))
    
    # This creates the Y-axis (counts)
    counts = zeros(length(edges) - 1)
    for t in times
        idx = floor(Int, (t - t_start) / binsize) + 1
        if 1 <= idx <= length(counts)
            counts[idx] += 1
        end
    end

    # Return the X (bin starts) and Y (counts)
    return edges[1:end-1], counts
end

# This tells Julia how to plot a LightCurve object
@recipe function f(lc::LightCurve; show_errors=true)
    title  --> "LightCurve (bin: $(lc.dt)s)"
    xlabel --> "Time (s)"
    ylabel --> "Counts"
    seriestype --> :steppre

    # FIX: Use lc.count_error instead of lc.errors
    if show_errors && hasfield(typeof(lc), :count_error) && !isnothing(lc.count_error)
        yerror := lc.count_error
    end

    # Handle GTI shading for LightCurves
    # In your version, we check if GTI exists in the extra metadata
    if haskey(lc.metadata.extra, "gti")
        gtis = lc.metadata.extra["gti"]
        for i in 1:size(gtis, 1)
            @series begin
                seriestype := :vspan
                fillcolor := :green
                fillalpha := 0.1
                label := ""
                primary := false
                [gtis[i, 1], gtis[i, 2]]
            end
        end
    end

    return lc.time, lc.counts
end
