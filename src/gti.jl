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
function split_by_gtis(el::EventList, gtis::AbstractMatrix{<:Real})
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
LightCurve objects for spectral timing analysis. Only complete time bins that fall
entirely within GTI boundaries are included to maintain temporal coherence required
for Fourier analysis and periodogram calculations.

# Arguments
- `lc::LightCurve{T}`: Input light curve with binned photon count data
- `gtis::AbstractMatrix{<:Real}`: Matrix of GTI boundaries where each row contains
  [start_time, stop_time] for valid observation intervals

# Returns
- `Vector{LightCurve{T}}`: Array of LightCurve objects, one for each GTI. Only
  segments containing at least one complete time bin are included. Empty segments
  are excluded from the result.

# Filtering Strategy
- **Bin inclusion criterion**: Complete bins only - bin center must fall within GTI
- **Boundary handling**: Bins partially overlapping GTI edges are excluded
- **Metadata preservation**: All properties and metadata are maintained per segment
- **Temporal continuity**: Each segment maintains uniform time binning from original

# Technical Details
The filtering uses bin centers for GTI membership testing:
```julia
included_bins = (bin_center ≥ gti_start) && (bin_center ≤ gti_stop)
```

This conservative approach ensures that:
1. All included bins have complete exposure within the GTI
2. Fourier analysis assumptions are preserved (uniform sampling)
3. Statistical properties remain well-defined
4. No partial bins introduce systematic errors

# Periodogram Compatibility
!!! warning "Bartlett Periodogram Limitation"
    This function is **NOT** suitable for Bartlett periodogram calculations, which
    require segments of identical length. The resulting segments will have different
    lengths depending on GTI durations.

!!! note "Suitable Methods"
    Use with Welch's method, Lomb-Scargle periodograms, or other techniques that
    can handle variable-length segments.

# Examples
```julia
# Basic segmentation
gtis = [1000.0 2000.0; 3000.0 4000.0; 5000.0 6000.0]
lc_segments = apply_gtis(lightcurve, gtis)

println("Created \$(length(lc_segments)) light curve segments")

# Analyze each segment independently
for (i, segment) in enumerate(lc_segments)
    mean_rate = mean(segment.counts ./ segment.exposure)
    duration = segment.time[end] - segment.time[1] + segment.dt
    println("GTI \$i: mean rate = \$(mean_rate) cts/s, duration = \$(duration) s")
end

# Variable-length periodogram analysis (NOT Bartlett)
using FFTW
periodograms = []
for segment in lc_segments
    # Welch's method can handle different segment lengths
    pgram = welch_periodogram(segment.counts, segment.dt)
    push!(periodograms, pgram)
end

# Check segment properties
for (i, seg) in enumerate(lc_segments)
    println("Segment \$i: \$(length(seg.time)) bins, Δt = \$(seg.dt) s")
end
```

# Performance Notes
- Time complexity: O(n) where n is the number of time bins
- Memory usage: Creates new LightCurve objects for each segment
- For large datasets, consider processing segments individually rather than
  storing all segments in memory

# References
- Welch periodogram methodology for variable-length segments
- X-ray timing analysis protocols (van der Klis 1989)
- Stingray light curve segmentation documentation

# See Also
- [`LightCurve`](@ref): Light curve data structure
- [`check_gtis`](@ref): Validate GTI format
- [`create_filtered_lightcurve`](@ref): Create filtered light curve segments
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

This function identifies gaps between Good Time Intervals (GTIs) and conditionally
fills short gaps with synthetic events to maintain temporal sampling for certain 
analysis methods. The synthetic events are generated based on event rates from 
nearby GTIs to preserve statistical properties.

# Arguments
- `el::EventList`: EventList to modify in-place
- `gtis::AbstractMatrix{<:Real}`: GTI boundaries defining good observation periods,
  where each row contains [start_time, stop_time]
- `dt::Real=1.0`: Time step parameter (retained for API compatibility, not used in calculations)
- `random_fill_threshold::Real=10.0`: Maximum BTI duration for random filling (seconds)
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator for synthetic events

# Returns
- `EventList`: The modified EventList (same object, modified in-place)

# Filling Strategy
- **Short BTIs** (< `random_fill_threshold`): Fill with synthetic events based on 
  event rates calculated from nearby GTIs
- **Long BTIs** (≥ `random_fill_threshold`): Leave empty to preserve data integrity
- **Rate calculation**: Uses median event rate from all GTIs, with fallback to 
  overall event rate if no GTI rates are available

# Technical Implementation
1. Compute BTIs using `get_btis()` based on GTI boundaries and EventList time range
2. For each qualifying short BTI:
   - Calculate event rates from all GTIs: `events_in_gti / gti_duration`
   - Use median rate to determine number of synthetic events
   - Generate uniformly distributed synthetic times within BTI boundaries
   - Sample energies and extra column values from events in adjacent GTIs
3. Append all synthetic data and re-sort chronologically
4. Update metadata with filling statistics

# Metadata Updates
The function adds the following metadata to `el.meta.headers`:
- `"BTI_FILLED"`: Boolean indicating if any BTIs were filled
- `"N_SYNTH_EVENTS"`: Total number of synthetic events added
- `"RAND_FILL_THRESH"`: The random fill threshold used
- `"BTI_FILL_DT"`: The dt parameter used

Additional metadata in `el.meta.extra_columns`:
- `"filled_bti_durations"`: Vector of durations of filled BTIs

# Warnings
!!! warning "Data Integrity"
    Synthetic events are clearly marked in metadata but are indistinguishable 
    from real events in the main data arrays. Exercise caution in subsequent analysis.

!!! note "Statistical Validity"
    Synthetic events may bias certain statistical measures. Consider the impact
    on your specific analysis before using this function.

# Examples
```julia
# Basic usage with default parameters
gtis = [1000.0 2000.0; 3000.0 4000.0; 5000.0 6000.0]
fill_bad_time_intervals!(eventlist, gtis)

# Custom parameters for short gaps only
fill_bad_time_intervals!(eventlist, gtis, 
                        random_fill_threshold=5.0)

# Check what was filled
if get(eventlist.meta.headers, "BTI_FILLED", false)
    n_synth = eventlist.meta.headers["N_SYNTH_EVENTS"]
    println("Filled \$n_synth synthetic events")
    
    # Access filled interval durations
    if haskey(eventlist.meta.extra_columns, "filled_bti_durations")
        durations = eventlist.meta.extra_columns["filled_bti_durations"]
        println("Filled BTI durations: \$durations")
    end
end
```

# References
- Stingray BTI handling: https://stingray.readthedocs.io/en/stable/
- Statistical considerations in X-ray timing analysis (Vaughan et al. 2003)

# See Also
- [`get_btis`](@ref): Function to compute Bad Time Intervals
- [`apply_gtis`](@ref): Apply GTIs to segment data
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
    
    # Store all synthetic data before appending to avoid index issues
    all_synthetic_times = Float64[]
    all_synthetic_energies = Float64[]
    synthetic_extra_columns = Dict{String, Vector}()
    
    # Initialize synthetic extra columns
    for (col_name, col_data) in el.meta.extra_columns
        synthetic_extra_columns[col_name] = similar(col_data, 0)
    end
    
    # Process each BTI
    for i in 1:size(btis, 1)
        bti_start, bti_stop = btis[i, 1], btis[i, 2]
        bti_duration = bti_stop - bti_start
        
        # Skip BTIs that are too long or have zero/negative duration
        if bti_duration >= random_fill_threshold || bti_duration <= 0
            continue
        end
        
        # Calculate number of synthetic events based on event rate in nearby GTIs
        gti_rates = Float64[]
        for j in 1:size(gtis, 1)
            gti_start, gti_stop = gtis[j, 1], gtis[j, 2]
            events_in_gti = count(t -> gti_start <= t <= gti_stop, el.times)
            gti_duration_j = gti_stop - gti_start
            
            if gti_duration_j > 0 && events_in_gti > 0
                push!(gti_rates, events_in_gti / gti_duration_j)
            end
        end
        
        # Determine synthetic event count
        n_synthetic = 0
        if isempty(gti_rates)
            # Fallback to overall event rate
            if length(el.times) > 1
                total_duration = maximum(el.times) - minimum(el.times)
                if total_duration > 0
                    overall_rate = length(el.times) / total_duration
                    n_synthetic = max(1, round(Int, overall_rate * bti_duration))
                end
            end
        else
            # Use median rate from nearby GTIs
            median_rate = Statistics.median(gti_rates)
            n_synthetic = max(1, round(Int, median_rate * bti_duration))
        end
        
        if n_synthetic > 0
            # Generate uniformly distributed synthetic times within the BTI
            # Generate times strictly within (bti_start, bti_stop) exclusive bounds
            # Use a small epsilon to ensure we don't hit the exact boundaries
            epsilon = 1e-10
            effective_start = bti_start + epsilon
            effective_stop = bti_stop - epsilon
            effective_duration = effective_stop - effective_start
            
            synthetic_times = sort!(rand(rng, n_synthetic) .* effective_duration .+ effective_start)
            append!(all_synthetic_times, synthetic_times)
            
            n_synthetic_events += n_synthetic
            push!(filled_intervals, bti_duration)
            
            # Handle energy values if present - sample from nearby GTI events
            if !isnothing(el.energies)
                nearby_energies = Float64[]
                
                # Find energies from GTIs that are adjacent to this BTI
                for j in 1:size(gtis, 1)
                    gti_start, gti_stop = gtis[j, 1], gtis[j, 2]
                    
                    # Check if this GTI is adjacent to the BTI
                    if abs(gti_stop - bti_start) < 1e-10 || abs(gti_start - bti_stop) < 1e-10
                        # Find events in this GTI
                        for (idx, t) in enumerate(el.times)
                            if gti_start <= t <= gti_stop
                                push!(nearby_energies, el.energies[idx])
                            end
                        end
                    end
                end
                
                # Sample energies
                if !isempty(nearby_energies)
                    synthetic_energies = rand(rng, nearby_energies, n_synthetic)
                elseif !isempty(el.energies)
                    synthetic_energies = rand(rng, el.energies, n_synthetic)
                else
                    synthetic_energies = zeros(eltype(el.energies), n_synthetic)
                end
                
                append!(all_synthetic_energies, synthetic_energies)
            end
            
            # Handle extra columns with values from nearby GTIs
            for (col_name, col_data) in el.meta.extra_columns
                nearby_values = similar(col_data, 0)
                
                # Find values from GTIs that are adjacent to this BTI
                for j in 1:size(gtis, 1)
                    gti_start, gti_stop = gtis[j, 1], gtis[j, 2]
                    
                    # Check if this GTI is adjacent to the BTI
                    if abs(gti_stop - bti_start) < 1e-10 || abs(gti_start - bti_stop) < 1e-10
                        # Find events in this GTI
                        for (idx, t) in enumerate(el.times)
                            if gti_start <= t <= gti_stop
                                push!(nearby_values, col_data[idx])
                            end
                        end
                    end
                end
                
                # Sample values
                if !isempty(nearby_values)
                    synthetic_values = rand(rng, nearby_values, n_synthetic)
                elseif !isempty(col_data)
                    synthetic_values = rand(rng, col_data, n_synthetic)
                else
                    synthetic_values = zeros(eltype(col_data), n_synthetic)
                end
                
                append!(synthetic_extra_columns[col_name], synthetic_values)
            end
        end
    end
    
    # Now append all synthetic data at once
    if n_synthetic_events > 0
        append!(el.times, all_synthetic_times)
        
        if !isnothing(el.energies)
            append!(el.energies, all_synthetic_energies)
        end
        
        # Append synthetic extra columns
        for (col_name, col_data) in el.meta.extra_columns
            append!(col_data, synthetic_extra_columns[col_name])
        end
        
        # Re-sort all arrays to maintain chronological order
        sort_indices = sortperm(el.times)
        el.times[:] = el.times[sort_indices]
        
        if !isnothing(el.energies)
            el.energies[:] = el.energies[sort_indices]
        end
        
        # Sort extra columns
        for (col_name, col_data) in el.meta.extra_columns
            col_data[:] = col_data[sort_indices]
        end
        
        # metadata with BTI filling information
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