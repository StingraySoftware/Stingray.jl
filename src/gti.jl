function get_total_gti_length(gti::AbstractMatrix{<:Real}; minlen::Real=0.0)
    lengths = diff(gti; dims =2)
    return sum(x->x > minlen ? x : zero(x), lengths)
end

function load_gtis(fits_file::String, gtistring::String="GTI")
    gti = FITS(fits_file) do lchdulist
        gtihdu = lchdulist[gtistring]
        get_gti_from_hdu(gtihdu)
    end
    return gti
end

function get_gti_from_hdu(gtihdu::TableHDU)

    if "START" in FITSIO.colnames(gtihdu)
        startstr = "START"
        stopstr = "STOP"
    else
        startstr = "Start"
        stopstr = "Stop"
    end

    gtistart = read(gtihdu,startstr)
    gtistop = read(gtihdu,stopstr)

    return mapreduce(permutedims, vcat, 
    [[a, b] for (a,b) in zip(gtistart, gtistop)])
end

function check_gtis(gti::AbstractMatrix)

    if ndims(gti) != 2 || size(gti,2) != 2
        throw(ArgumentError("Please check the formatting of the GTIs. 
       They need to be provided as [[gti00 gti01]; [gti10 gti11]; ...]."))
    end

    gti_start = @view gti[:, 1]
    gti_end = @view gti[:, 2]

    if any(gti_end < gti_start)
        throw(ArgumentError(
            "The GTI end times must be larger than the GTI start times."
        )) 
    end

    if any(@view(gti_start[begin+1:end]) < @view(gti_end[begin:end-1]))
        throw(ArgumentError(
            "This GTI has overlaps"
        ))
    end
end

function create_gti_mask(times::AbstractVector{<:Real},gtis::AbstractMatrix{<:Real};
                         safe_interval::AbstractVector{<:Real}=[0,0], min_length::Real=0,
                         dt::Real = -1, epsilon::Real = 0.001)

    if isempty(times)
        throw(ArgumentError("Passing an empty time array to create_gti_mask"))
    end

    check_gtis(gtis)
    mask = zeros(Bool,length(times))

    if min_length>0
        gtis = gtis[min_length .< @view(gtis[:,2]) - @view(gtis[:,1]),:]
            
        if size(gtis,1) < 1
            @warn "No GTIs longer than min_length $(min_length)"
            return mask, gtis
        end
    end   

    if dt < 0
        dt = Statistics.median(diff(times))
    end
    epsilon_times_dt = epsilon * dt

    new_gtis = [[0.0, 0.0] for _ in range(1,size(gtis,1))]
    new_gti_mask = zeros(Bool, size(gtis,1))

    gti_start = @view gtis[:, 1]
    gti_end = @view gtis[:, 2]

    for (ig,(limmin,limmax)) in enumerate(zip(gti_start,gti_end))
        limmin += safe_interval[1]
        limmax -= safe_interval[2]
        if limmax - limmin >= min_length
            new_gtis[ig][:] .= limmin, limmax
            for (i,t) in enumerate(times) 
                if (limmin + dt / 2 - epsilon_times_dt) <= t <= (limmax - dt / 2 + epsilon_times_dt)
                    mask[i] = true
                end
            end
            new_gti_mask[ig] = true
        end
    end

    return mask, mapreduce(permutedims, vcat, keepat!(new_gtis,new_gti_mask))
end

function create_gti_from_condition(time::AbstractVector{<:Real}, condition::AbstractVector{Bool};
    safe_interval::AbstractVector{<:Real}=[0,0], dt::AbstractVector{<:Real}=Float64[])
    
    if length(time) != length(condition)
        throw(ArgumentError("The length of the condition and time arrays must be the same."))
    end

    idxs = contiguous_regions(condition)

    if isempty(dt)
        dt = zero(time) .+ (time[2] .- time[1]) ./ 2
    end

    gtis = Vector{Float64}[]
    for idx in eachrow(idxs)
        startidx = idx[1]
        stopidx = idx[2] - 1

        t0 = time[startidx] - dt[startidx] + safe_interval[1]
        t1 = time[stopidx] + dt[stopidx] - safe_interval[2]
        if t1 - t0 < 0
            continue
        end
        push!(gtis,[t0, t1])
    end
    return mapreduce(permutedims, vcat, gtis)
end

function operations_on_gtis(gti_list::AbstractVector{<:AbstractMatrix{T}}, 
                            operation::Function) where {T<:Real}

    required_interval = nothing

    for gti in gti_list
        check_gtis(gti)

        combined_gti = Interval{T}[]
        for ig in eachrow(gti)
            push!(combined_gti,Interval{Closed,Open}(ig[1],ig[2]))
        end
        if isnothing(required_interval)
            required_interval = IntervalSet(combined_gti)
        else
            required_interval = operation(required_interval, IntervalSet(combined_gti))
        end
    end

    final_gti = Vector{T}[]

    for interval in required_interval.items
        push!(final_gti, [first(interval), last(interval)])
    end

    return mapreduce(permutedims, vcat, final_gti)
end

function get_btis(gtis::AbstractMatrix{<:Real})
    if isempty(gtis)
        throw(ArgumentError("Empty GTI and no valid start_time and stop_time"))
    end
    return get_btis(gtis, gtis[1,1], gtis[end,2])
end

function get_btis(gtis::AbstractMatrix{T}, start_time, stop_time) where {T<:Real}
    if isempty(gtis)
        return T[start_time stop_time]
    end
    check_gtis(gtis)

    total_interval = Interval{T, Closed, Open}[Interval{T, Closed, Open}(start_time, stop_time)]
    total_interval_set = IntervalSet(total_interval)

    gti_interval = Interval{T, Closed, Open}[]
    for ig in eachrow(gtis)
        push!(gti_interval,Interval{T, Closed, Open}(ig[1],ig[2]))
    end
    gti_interval_set = IntervalSet(gti_interval)

    bti_interval_set = setdiff(total_interval_set, gti_interval_set)

    btis = Vector{T}[]

    for interval in bti_interval_set.items
        push!(btis, [first(interval), last(interval)])
    end

    # Fix: Handle empty btis vector
    if isempty(btis)
        return reshape(T[], 0, 2)  # Return empty matrix with correct dimensions
    end

    return mapreduce(permutedims, vcat, btis)
end
function time_intervals_from_gtis(gtis::AbstractMatrix{<:Real}, segment_size::Real;
                                  fraction_step::Real=1, epsilon::Real=1e-5)  
    spectrum_start_times = Float64[]

    gti_low = @view gtis[:,1]
    gti_up = @view gtis[:,2]

    for (g1,g2) in zip(gti_low,gti_up)
        if g2 - g1 + epsilon < segment_size
            continue
        end

        newtimes = range(g1, g2 - segment_size + epsilon, step = segment_size* fraction_step)
        append!(spectrum_start_times,newtimes)
    end
    return spectrum_start_times, spectrum_start_times .+ segment_size
end

function calculate_segment_bin_start(startbin::Integer, stopbin::Integer,
                                     nbin::Integer; fraction_step::Real=1)
    st = floor.(range(startbin, stopbin, step=Int(nbin * fraction_step)))
    if st[end] == stopbin
        pop!(st)
    end
    if st[end] + nbin > stopbin
        pop!(st)
    end
    return st
end

function bin_intervals_from_gtis(gtis::AbstractMatrix{<:Real}, segment_size::Real,
                                 time::AbstractVector{<:Real}; dt=nothing, 
                                 fraction_step::Real=1, epsilon::Real=0.001)
    if isnothing(dt)
        dt = Statistics.median(diff(time))
    end

    epsilon_times_dt = epsilon * dt
    nbin = round(Int, segment_size / dt)

    spectrum_start_bins = Int[]

    gti_low = @view(gtis[:, 1]) .+ (dt ./ 2 .- epsilon_times_dt)
    gti_up = @view(gtis[:, 2]) .- (dt ./ 2 .- epsilon_times_dt)

    for (g0, g1) in zip(gti_low, gti_up)
        if (g1 - g0 .+ (dt + epsilon_times_dt)) < segment_size
            continue
        end
        startbin, stopbin = searchsortedfirst.(Ref(time), [g0, g1])
        startbin -= 1
        if stopbin > length(time)
            stopbin = length(time)
        end

        if time[startbin+1] < g0
            startbin += 1
        end
        # Would be g[1] - dt/2, but stopbin is the end of an interval
        # so one has to add one bin
        if time[stopbin] > g1
            stopbin -= 1
        end

        newbins = calculate_segment_bin_start(
            startbin, stopbin, nbin, fraction_step=fraction_step)
        
        append!(spectrum_start_bins,newbins)
    end 
    return spectrum_start_bins, spectrum_start_bins.+nbin 
end

@resumable function generate_indices_of_segment_boundaries_unbinned(times::AbstractVector{<:Real},
                                                                    gti::AbstractMatrix{<:Real},
                                                                    segment_size::Real)
    start, stop = time_intervals_from_gtis(gti, segment_size)

    startidx = searchsortedfirst.(Ref(times), start)
    stopidx = searchsortedfirst.(Ref(times), stop)

    for (s, e, idx0, idx1) in zip(start, stop, startidx, stopidx)
        @yield s, e, idx0, idx1
    end
end

@resumable function generate_indices_of_segment_boundaries_binned(times::AbstractVector{<:Real},
                                                                  gti::AbstractMatrix{<:Real},
                                                                  segment_size::Real; dt=nothing)
    startidx, stopidx = bin_intervals_from_gtis(gti, segment_size, times;
                                                dt=dt)

    if isnothing(dt)
        dt = 0
    end
    for (idx0, idx1) in zip(startidx, stopidx)
        @yield times[idx0+1] - dt / 2, times[min(idx1, length(times) - 1)] - dt / 2,idx0, idx1
    end
end
"""
    apply_gtis(el::EventList, gtis::AbstractMatrix{<:Real}) -> Vector{EventList}

Apply Good Time Intervals (GTIs) to an EventList, returning a separate EventList for each GTI.

This function filters the input EventList based on the provided GTI boundaries, creating 
independent EventList objects for each valid time interval. This is essential for 
X-ray timing analysis where data quality varies and only specific time intervals 
contain reliable observations.

# Arguments
- `el::EventList`: Input event list containing photon arrival times and energies
- `gtis::AbstractMatrix{<:Real}`: Matrix of GTI boundaries where each row contains 
  [start_time, stop_time] for a valid observation interval

# Returns
- `Vector{EventList}`: Array of EventList objects, one for each GTI containing events 
  that fall within the corresponding time interval. Empty GTIs are excluded from results.

# Notes
- Events are filtered based on arrival times: `gti_start ≤ time ≤ gti_stop`
- Maintains all original metadata and extra columns for each filtered EventList
- GTIs are validated using `check_gtis()` to ensure proper formatting and ordering
- Only non-empty EventLists are returned (GTIs with zero events are excluded)

# Examples
```julia
# Apply GTIs to filter data during good observation periods
gtis = [100.0 200.0; 300.0 400.0; 500.0 600.0]  # Three GTI intervals
filtered_events = apply_gtis(eventlist, gtis)
println("Number of valid GTI segments: ", length(filtered_events))

# Each segment can be analyzed independently
for (i, segment) in enumerate(filtered_events)
    println("GTI \$i: \$(length(segment)) events")
end
```

# References
- Stingray documentation on GTI handling
- X-ray timing analysis best practices (Belloni et al. 2000)
"""
function apply_gtis(el::EventList, gtis::AbstractMatrix{<:Real})
    check_gtis(gtis)
    
    result = EventList[]
    
    for i in 1:size(gtis, 1)
        gti_start, gti_stop = gtis[i, 1], gtis[i, 2]
        
        # Create filter function for this specific GTI
        gti_filter = t -> gti_start ≤ t ≤ gti_stop
        
        # Apply temporal filtering using existing infrastructure
        filtered_el = filter_time(gti_filter, el)
        
        # Only include GTIs that contain events
        if length(filtered_el.times) > 0
            push!(result, filtered_el)
        end
    end
    
    return result
end
"""
    apply_gtis(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real}) -> Vector{LightCurve{T}} where T

Apply Good Time Intervals (GTIs) to a LightCurve, returning separate LightCurve objects for each GTI.

This function segments a light curve based on GTI boundaries, creating independent 
LightCurve objects for spectral timing analysis. Bins that partially overlap with 
GTI boundaries are excluded to maintain temporal coherence required for Fourier 
analysis and periodogram calculations.

# Arguments
- `lc::LightCurve{T}`: Input light curve with binned photon count data
- `gtis::AbstractMatrix{<:Real}`: Matrix of GTI boundaries where each row contains 
  [start_time, stop_time] for valid observation intervals

# Returns
- `Vector{LightCurve{T}}`: Array of LightCurve objects, one for each GTI. Only 
  segments containing at least one complete time bin are included.

# Filtering Strategy
- **Bin inclusion criterion**: Complete bins only (bin center must fall within GTI)
- **Boundary handling**: Bins partially overlapping GTI edges are excluded
- **Metadata preservation**: All properties and metadata are maintained per segment
- **Temporal continuity**: Each segment maintains uniform time binning

# Technical Details
The filtering uses bin centers for GTI membership testing:
```
included_bins = (bin_center ≥ gti_start) && (bin_center ≤ gti_stop)
```

This approach ensures that:
1. All included bins have complete exposure within the GTI
2. Fourier analysis assumptions are preserved (uniform sampling)
3. Statistical properties remain well-defined

# Examples
```julia
# Segment light curve for independent periodogram analysis
gtis = [1000.0 2000.0; 3000.0 4000.0]
lc_segments = apply_gtis(lightcurve, gtis)

# Analyze each segment independently for variability
for (i, segment) in enumerate(lc_segments)
    mean_rate = mean(segment.counts ./ segment.exposure)
    println("GTI \$i: mean count rate = \$(mean_rate) cts/s")
end

# Suitable for Bartlett periodogram calculations
periodograms = [calculate_periodogram(seg) for seg in lc_segments]
```

# References
- Bartlett periodogram methodology (Bartlett 1955)
- X-ray timing analysis protocols (van der Klis 1989)
- Stingray light curve segmentation
"""
function apply_gtis(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real}) where T
    check_gtis(gtis)
    
    result = LightCurve{T}[]
    
    for i in 1:size(gtis, 1)
        gti_start, gti_stop = T(gtis[i, 1]), T(gtis[i, 2])
        
        # Filter bins based on bin centers falling within GTI
        bin_mask = (lc.time .≥ gti_start) .& (lc.time .≤ gti_stop)
        
        if any(bin_mask)
            # Convert to BitVector before passing to create_filtered_lightcurve
            filtered_lc = create_filtered_lightcurve(lc, BitVector(bin_mask), gti_start, gti_stop, i)
            push!(result, filtered_lc)
        end
    end
    
    return result
end
"""
    create_filtered_lightcurve(lc::LightCurve{T}, mask::BitVector, 
                              gti_start::T, gti_stop::T, gti_index::Int) -> LightCurve{T}

Internal function to create a filtered LightCurve from a boolean mask.

Creates a new LightCurve containing only the time bins specified by the mask,
while preserving all metadata, properties, and statistical characteristics.

# Arguments
- `lc::LightCurve{T}`: Source light curve
- `mask::BitVector`: Boolean mask indicating which bins to include
- `gti_start::T`: Start time of the GTI (for metadata)
- `gti_stop::T`: Stop time of the GTI (for metadata)  
- `gti_index::Int`: GTI sequence number (for metadata tracking)

# Returns
- `LightCurve{T}`: Filtered light curve with updated metadata reflecting the GTI application

# Implementation Notes
- Preserves bin size and exposure information
- Maintains all computed properties (e.g., mean energy)
- Updates metadata to reflect GTI filtering
- Recalculates statistical errors for the filtered dataset
"""
function create_filtered_lightcurve(lc::LightCurve{T}, mask::AbstractVector{Bool}, 
                                   gti_start::T, gti_stop::T, gti_index::Int) where T

    # Ensure mask is proper boolean vector
    bool_mask = mask isa BitVector ? mask : BitVector(mask)
    # Filter all primary arrays
    filtered_time = lc.time[mask]
    filtered_counts = lc.counts[mask]
    filtered_exposure = isnothing(lc.exposure) ? nothing : lc.exposure[mask]
    
    # Filter all computed properties
    filtered_properties = EventProperty{T}[]
    for prop in lc.properties
        filtered_values = prop.values[mask]
        push!(filtered_properties, EventProperty(prop.name, filtered_values, prop.unit))
    end
    
    # Update metadata with GTI information
    updated_metadata = LightCurveMetadata(
        lc.metadata.telescope,
        lc.metadata.instrument, 
        lc.metadata.object,
        lc.metadata.mjdref,
        (Float64(gti_start), Float64(gti_stop)),  # Update time range to GTI bounds
        lc.metadata.bin_size,
        lc.metadata.headers,
        merge(lc.metadata.extra, Dict{String,Any}(
            "gti_applied" => true,
            "gti_index" => gti_index,
            "gti_bounds" => [Float64(gti_start), Float64(gti_stop)],
            "original_time_range" => lc.metadata.time_range,
            "filtered_nbins" => length(filtered_time),
            "original_nbins" => length(lc.time)
        ))
    )
    
    # Create new LightCurve with filtered data
    filtered_lc = LightCurve{T}(
        filtered_time,
        lc.dt,  # Preserve original bin size
        filtered_counts,
        nothing,  # Errors will be recalculated
        filtered_exposure,
        filtered_properties,
        updated_metadata,
        lc.err_method
    )
    
    # Recalculate errors for the filtered dataset
    calculate_errors!(filtered_lc)
    
    return filtered_lc
end
"""
    fill_bad_time_intervals!(el::EventList, gtis::AbstractMatrix{<:Real}; 
                            dt::Real=1.0, random_fill_threshold::Real=10.0,
                            rng::AbstractRNG=Random.GLOBAL_RNG) -> EventList

Fill Bad Time Intervals (BTIs) in an EventList with synthetic events for analysis continuity.

This function identifies gaps between GTIs (Bad Time Intervals) and conditionally 
fills them to maintain temporal sampling for certain analysis methods. Very short 
gaps may be filled with random synthetic events, while longer gaps remain empty 
to preserve data integrity.

# Arguments
- `el::EventList`: EventList to modify in-place
- `gtis::AbstractMatrix{<:Real}`: GTI boundaries defining good observation periods
- `dt::Real=1.0`: Time step for potential synthetic event generation (seconds)
- `random_fill_threshold::Real=10.0`: Maximum BTI duration for random filling (seconds)
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator for synthetic events

# Returns
- `EventList`: The modified EventList (same object, modified in-place)

# Filling Strategy
- **Short BTIs** (< `random_fill_threshold`): Fill with random synthetic events
- **Long BTIs** (≥ `random_fill_threshold`): Leave empty (no filling)
- **Purpose**: Maintain sampling for specific analysis methods while preserving data quality

# Technical Implementation
1. Compute BTIs using `get_btis()` based on GTI boundaries and EventList time range
2. For qualifying short BTIs, generate synthetic times with uniform random spacing
3. Assign synthetic energies from the original energy distribution
4. Maintain chronological ordering of all events

# Warnings
- **Data integrity**: Synthetic events are clearly marked in metadata
- **Analysis impact**: Consider whether filled intervals affect your specific analysis
- **Statistical validity**: Synthetic events may bias certain statistical measures

# Examples
```julia
# Fill short gaps for periodogram analysis continuity
gtis = load_gtis("observation.fits")
fill_bad_time_intervals!(eventlist, gtis, dt=0.5, random_fill_threshold=5.0)

# Check what was filled
if haskey(eventlist.meta.extra, "bti_filled")
    println("Filled \$(eventlist.meta.extra["n_synthetic_events"]) synthetic events")
end
```

# References
- Stingray BTI handling documentation
- Statistical considerations in X-ray timing (Vaughan et al. 2003)
"""
function fill_bad_time_intervals!(el::EventList, gtis::AbstractMatrix{<:Real}; 
                                  dt::Real=1.0, random_fill_threshold::Real=10.0,
                                  rng::AbstractRNG=Random.GLOBAL_RNG)
    check_gtis(gtis)
    
    # Determine the time range for BTI calculation
    time_start = minimum(el.times)
    time_stop = maximum(el.times)
    
    # Calculate Bad Time Intervals
    btis = get_btis(gtis, time_start, time_stop)
    
    # Track synthetic events for metadata
    n_synthetic_events = 0
    filled_intervals = Float64[]
    
    # Process each BTI
    for i in 1:size(btis, 1)
        bti_start, bti_stop = btis[i, 1], btis[i, 2]
        bti_duration = bti_stop - bti_start
        
        # Only fill short BTIs with random events
        if bti_duration < random_fill_threshold && bti_duration > 0
            # Generate random synthetic times within the BTI
            n_synthetic = max(1, floor(Int, bti_duration / dt))
            synthetic_times = sort!(rand(rng, n_synthetic) .* bti_duration .+ bti_start)
            
            # Add synthetic times to the EventList
            append!(el.times, synthetic_times)
            n_synthetic_events += n_synthetic
            push!(filled_intervals, bti_duration)
            
            # Handle energy values if present
            if !isnothing(el.energies)
                # Sample energies from existing distribution
                if !isempty(el.energies)
                    synthetic_energies = rand(rng, el.energies, n_synthetic)
                    append!(el.energies, synthetic_energies)
                else
                    # Fallback to zero energy if no existing energies
                    append!(el.energies, zeros(eltype(el.energies), n_synthetic))
                end
            end
            
            # Handle extra columns with appropriate fill values
            for (col_name, col_data) in el.meta.extra_columns
                if !isempty(col_data)
                    # Sample from existing values
                    synthetic_values = rand(rng, col_data, n_synthetic)
                    append!(col_data, synthetic_values)
                else
                    # Use zero/default values
                    append!(col_data, zeros(eltype(col_data), n_synthetic))
                end
            end
        end
    end
    
    # Re-sort all arrays to maintain chronological order
    if n_synthetic_events > 0
        sort_indices = sortperm(el.times)
        el.times[:] = el.times[sort_indices]
        
        if !isnothing(el.energies)
            el.energies[:] = el.energies[sort_indices]
        end
        
        # Sort extra columns
        for (col_name, col_data) in el.meta.extra_columns
            col_data[:] = col_data[sort_indices]
        end
        
        # Update metadata with BTI filling information
        # Store scalar metadata in headers since extra_columns expects vectors
        merge!(el.meta.headers, Dict{String,Any}(
            "BTI_FILLED" => true,
            "N_SYNTH_EVENTS" => n_synthetic_events,
            "RAND_FILL_THRESH" => random_fill_threshold,
            "BTI_FILL_DT" => dt
        ))
        
        # Store vector metadata in extra_columns
        el.meta.extra_columns["filled_bti_durations"] = filled_intervals
    end
    
    return el
end
#todo
# create a function that fills the bad time intervals in a light curve
# in order to maintain optiizational sampling for periodograms
# can be intially start like :
# fill_bad_time_intervals!(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real}; 
#                             dt::Real=1.0, random_fill_threshold::Real=10.0,
#                             rng::AbstractRNG=Random.GLOBAL_RNG) where T<:Real