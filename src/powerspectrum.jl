"""
    PowerSpectrum{T}

$(TYPEDEF)

Power spectrum for a single light curve segment.

$(TYPEDFIELDS)

# Examples
```julia
# Create power spectrum from light curve
lc = LightCurve(times, counts)
ps = Powerspectrum(lc, norm="leahy")

# Access properties
println(ps.freq)    # Frequency array
println(ps.power)   # Power values
println(ps.norm)    # Normalization type
```
"""
struct PowerSpectrum{T}
    "Frequencies in Hz"
    freq::Vector{T}
    "Power values in requested normalization"
    power::Vector{T}
    "Power errors based on normalization type"
    power_err::Vector{T}
    "Normalization type (leahy, frac, rms, abs)"
    norm::String
    "Frequency resolution (Hz)"
    df::T
    "Total number of photons"
    nphots::Int
    "Number of segments (1 for single spectrum)"
    m::Int
    "Number of frequencies"
    n::Int
    "Original light curve metadata"
    metadata::Union{LightCurveMetadata, FITSMetadata}
end

"""
    AveragedPowerspectrum{T}

$(TYPEDEF)

Averaged power spectrum from multiple light curve segments.

$(TYPEDFIELDS)

# Examples
```julia
# Create averaged power spectrum from light curve
lc = LightCurve(times, counts)
ps_avg = AveragedPowerspectrum(lc, 1024.0, norm="leahy")

# Access properties
println(ps_avg.freq)         # Frequency array
println(ps_avg.power)        # Averaged power values
println(ps_avg.segment_size) # Size of segments used
```
"""
struct AveragedPowerspectrum{T}
    "Frequencies in Hz"
    freq::Vector{T}
    "Averaged power values"
    power::Vector{T}
    "Errors on averaged powers"
    power_err::Vector{T}
    "Normalization type (leahy, frac, rms, abs)"
    norm::String
    "Frequency resolution (Hz)"
    df::T
    "Size of segments in seconds"
    segment_size::T
    "Total number of photons"
    nphots::Int
    "Number of segments averaged"
    m::Int
    "Mean count rate (cts/s)"
    mean_rate::T
    "Number of frequencies"
    n::Int
    "Original light curve metadata"
    metadata::Union{LightCurveMetadata, FITSMetadata}
end

#lightcurve==>

"""
    Powerspectrum(lc::LightCurve{T}; norm::String="frac") where T<:Real

Create power spectrum from a light curve.

# Arguments
- `lc`: Light curve to analyze
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")

# Returns
- `PowerSpectrum` object

# Examples
```julia
ps = Powerspectrum(lc, norm="leahy")
```
"""
function Powerspectrum(lc::LightCurve{T}; norm::String="frac") where T<:Real
    bin_size = lc.metadata.bin_size
    n_bins = length(lc.counts)
    
    # Input validation
    n_bins > 1 || throw(ArgumentError("Light curve must have more than 1 bin"))
    bin_size > 0 || throw(ArgumentError("Bin size must be positive"))
    
    # Calculate FFT
    ft = fft(lc.counts)
    
    # Get frequency array in Hz - use proper sampling frequency
    freqs = fftfreq(n_bins, 1/bin_size)
    pos_freq_idx = positive_fft_bins(n_bins)
    freqs = freqs[pos_freq_idx]
    
    # Calculate power
    unnorm_power = abs2.(ft[pos_freq_idx])
    
    # Normalize power
    power = normalize_periodograms(
        unnorm_power,
        bin_size,
        n_bins;
        mean_flux = mean(lc.counts),
        n_ph = sum(lc.counts),
        norm = norm,
        power_type = "all"
    )
    
    # Calculate errors based on normalization
    power_err = if norm == "leahy"
        fill(2.0, length(power))
    elseif norm in ["frac", "rms"]
        power ./ sqrt(1)
    else
        sqrt.(power)
    end
    
    return PowerSpectrum{T}(
        freqs,
        power,
        power_err,
        norm,
        freqs[2] - freqs[1],
        sum(lc.counts),
        1,
        length(freqs),
        lc.metadata
    )
end

"""
    AveragedPowerspectrum(lc::LightCurve{T}, segment_size::Real; norm::String="frac", epsilon::Real=1e-5) where T<:Real

Create averaged power spectrum from a light curve divided into segments.

# Arguments
- `lc`: Light curve to analyze
- `segment_size`: Size of segments in seconds
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")
- `epsilon`: Tolerance for segment boundaries

# Returns
- `AveragedPowerspectrum` object

# Examples
```julia
ps_avg = AveragedPowerspectrum(lc, 1024.0, norm="leahy")
```
"""
function AveragedPowerspectrum(lc::LightCurve{T}, segment_size::Real; 
                             norm::String="frac", 
                             epsilon::Real=1e-5) where T<:Real
    # this segment is for debug/valid can be removed build for [Comprehensive Error Handling tests]
    if isnan(segment_size)
        throw(ArgumentError("Segment size cannot be NaN"))
    end
    if isinf(segment_size)
        throw(ArgumentError("Segment size cannot be Inf"))
    end
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive"))
    end

    bin_size = lc.metadata.bin_size
    n_bins_per_segment = round(Int, segment_size / bin_size)

    # Get GTIs from metadata
    gtis = if hasfield(typeof(lc.metadata), :gti) && !isnothing(lc.metadata.gti)
        lc.metadata.gti
    elseif haskey(lc.metadata.extra, "gti")
        lc.metadata.extra["gti"]
    else
        throw(ArgumentError("No GTI information found in metadata"))
    end
    # UNBINNED segmentation like events method
    segment_generator = generate_indices_of_segment_boundaries_unbinned(
        lc.time, gtis, segment_size
    )
    
    # Initialize arrays
    freqs = fftfreq(n_bins_per_segment, 1/bin_size)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero=false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]
    
    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0
    
    # Process each segment
    for (start_time, stop_time, start_idx, stop_idx) in segment_generator
        
        # Find corresponding light curve bins
        bin_mask = (lc.time .>= start_time) .& (lc.time .< stop_time)
        segment_counts = lc.counts[bin_mask]
        
        actual_length = length(segment_counts)
        
        if actual_length < n_bins_per_segment * 0.8
            continue
        elseif actual_length > n_bins_per_segment
            # Truncate if too long
            segment_counts = segment_counts[1:n_bins_per_segment]
        elseif actual_length < n_bins_per_segment
            # Pad with zeros if too short (but still above threshold)
            padded_counts = zeros(Int, n_bins_per_segment)
            padded_counts[1:actual_length] .= segment_counts
            segment_counts = padded_counts
        end
        
        segment_sum = sum(segment_counts)
        if segment_sum == 0
            @debug "Skipping segment: zero counts"
            continue
        end
        
        # Calculate FFT and power
        ft = fft(segment_counts)
        unnorm_power = abs2.(ft[pos_freq_idx])
        
        # Normalize power
        power = normalize_periodograms(
            unnorm_power,
            bin_size,
            n_bins_per_segment;
            mean_flux = mean(segment_counts),
            n_ph = segment_sum,
            norm = norm,
            power_type = "all"
        )
        
        total_power .+= power
        total_counts += segment_sum
        n_segments_used += 1
        
        @debug "Processed segment" start_time stop_time n_counts=segment_sum mean_rate=segment_sum/segment_size
    end
    
    @debug "Segment processing complete" segments_used=n_segments_used
    
    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end
    
    # Calculate final results
    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)
    
    # Calculate errors
    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end
    
    return AveragedPowerspectrum{T}(
        freqs,
        avg_power,
        power_err,
        norm,
        df,
        segment_size,
        total_counts,
        n_segments_used,
        mean_rate,
        length(freqs),
        lc.metadata
    )
end

#Eventlist==>

"""
    Powerspectrum(events::EventList{Vector{T}, M}; norm::String="frac", dt::Real=1.0) where {T<:Real, M}

Create power spectrum from an event list by first binning the events.

# Arguments
- `events`: EventList to analyze
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")
- `dt`: Bin size in seconds for creating light curve

# Returns
- `PowerSpectrum` object

# Examples
```julia
events = readevents("data.fits")
ps = Powerspectrum(events, norm="leahy", dt=0.1)
```
"""
function Powerspectrum(events::EventList{Vector{T}, M}, 
                                dt::Real, 
                                segment_size::Real; 
                                norm::String="leahy") where {T<:Real, M}
    
    # Input validation
    length(events) > 1 || throw(ArgumentError("EventList must have more than 1 event"))
    dt > 0 || throw(ArgumentError("Bin size must be positive"))
    segment_size > dt || throw(ArgumentError("Segment size must be larger than bin size"))
    
    n_bins_per_segment = round(Int, segment_size / dt)
    
    # Get GTIs from events metadata
    gtis = if has_gti(events)
        events.meta.gti
    else
        # Create single GTI spanning the entire observation
        time_span = extrema(events.times)
        reshape([time_span[1], time_span[2]], 1, 2)
    end
    
    @debug "Direct events processing" total_events=length(events) segment_size=segment_size dt=dt bins_per_segment=n_bins_per_segment
    
    # Use unbinned segment generator - this preserves exact event timing
    segment_generator = generate_indices_of_segment_boundaries_unbinned(
        events.times, gtis, segment_size
    )
    
    # Initialize frequency arrays
    freqs = fftfreq(n_bins_per_segment, 1/dt)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero=false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]
    
    # Initialize accumulators
    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0
    
    # Process each segment
    for (start_time, stop_time, start_idx, stop_idx) in segment_generator
        
        # Extract events in this segment - preserve exact timing
        segment_event_times = if start_idx <= stop_idx && start_idx > 0 && stop_idx <= length(events.times)
            @view events.times[start_idx:stop_idx]
        else
            # Fallback to time-based filtering for edge cases
            filter(t -> start_time <= t < stop_time, events.times)
        end
        
        # Skip segments with too few events
        if length(segment_event_times) < 2
            @debug "Skipping segment: too few events" start_time stop_time n_events=length(segment_event_times)
            continue
        end
        
        # Create time grid for this segment
        time_grid = range(start_time, stop=stop_time, length=n_bins_per_segment+1)
        bin_centers = (time_grid[1:end-1] + time_grid[2:end]) / 2
        
        # Bin events directly without creating LightCurve object
        # This preserves more control over the binning process
        counts = zeros(Int, n_bins_per_segment)
        
        for event_time in segment_event_times
            # Find which bin this event belongs to
            bin_idx = searchsortedfirst(time_grid, event_time)
            if 1 <= bin_idx <= n_bins_per_segment
                counts[bin_idx] += 1
            end
        end
        
        segment_total_counts = sum(counts)
        
        # Skip if no counts in this segment
        if segment_total_counts == 0
            @debug "Skipping segment: zero counts after binning" start_time stop_time
            continue
        end
        
        # Calculate FFT directly on counts
        ft = fft(counts)
        unnorm_power = abs2.(ft[pos_freq_idx])
        
        # Normalize using your existing function
        power = normalize_periodograms(
            unnorm_power,
            dt,
            n_bins_per_segment;
            mean_flux = mean(counts),
            n_ph = segment_total_counts,
            norm = norm,
            power_type = "all"
        )
        
        # Accumulate results
        total_power .+= power
        total_counts += segment_total_counts
        n_segments_used += 1
        
        @debug "Processed segment directly" start_time stop_time n_events=length(segment_event_times) n_counts=segment_total_counts mean_rate=segment_total_counts/segment_size
    end
    
    @debug "Direct processing complete" segments_used=n_segments_used total_counts=total_counts
    
    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end
    
    # Calculate final averaged results
    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)
    
    # Calculate errors based on normalization
    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end
    
    # Create metadata preserving original event information
    result_metadata = create_powerspectrum_metadata(events, dt, segment_size)
    
    return AveragedPowerspectrum{T}(
        freqs,
        avg_power,
        power_err,
        norm,
        df,
        segment_size,
        total_counts,
        n_segments_used,
        mean_rate,
        length(freqs),
        result_metadata
    )
end

"""
    create_powerspectrum_metadata(events::EventList, dt::Real, segment_size::Real) -> LightCurveMetadata

Create metadata for power spectrum analysis results from EventList.

Extracts FITS header information and creates metadata structure documenting
the power spectrum analysis parameters and data provenance.

# Arguments
- `events::EventList`: Source event data with FITS headers
- `dt::Real`: Time bin size used in analysis (seconds)
- `segment_size::Real`: Segment duration for averaging (seconds)

# Returns
`LightCurveMetadata`: Metadata with analysis parameters and GTI information

# Examples
```julia
metadata = create_powerspectrum_metadata(events, 0.1, 1024.0)
"""
function create_powerspectrum_metadata(events::EventList, dt::Real, segment_size::Real)
    headers = events.meta.headers
    telescope = try
        get(headers, "TELESCOP", "")
    catch
        haskey(headers, "TELESCOP") ? headers["TELESCOP"] : ""
    end
    
    instrument = try
        get(headers, "INSTRUME", "")
    catch
        haskey(headers, "INSTRUME") ? headers["INSTRUME"] : ""
    end
    
    object_name = try
        get(headers, "OBJECT", "")
    catch
        haskey(headers, "OBJECT") ? headers["OBJECT"] : ""
    end
    
    mjdref = try
        get(headers, "MJDREF", 0.0)
    catch
        haskey(headers, "MJDREF") ? headers["MJDREF"] : 0.0
    end
    
    # Get GTIs
    gtis = has_gti(events) ? events.meta.gti : reshape([minimum(events.times), maximum(events.times)], 1, 2)
    
    return LightCurveMetadata(
        telescope,
        instrument,
        object_name,
        mjdref,
        (minimum(events.times), maximum(events.times)),
        Float64(dt),
        [Dict{String,Any}("TELESCOP" => telescope, "INSTRUME" => instrument, "OBJECT" => object_name, "MJDREF" => mjdref)],
        Dict(
            "analysis_method" => "direct_events_processing",
            "original_file" => events.meta.filepath,
            "original_hdu" => events.meta.hdu,
            "energy_units" => events.meta.energy_units,
            "n_original_events" => length(events),
            "segment_size" => segment_size,
            "time_resolution" => dt,
            "gti" => gtis,
            "original_fits_header" => headers
        )
    )
end
"""
    AveragedPowerspectrum(events::EventList{Vector{T}, M}, segment_size::Real; 
                         norm::String="frac", dt::Real=1.0, 
                         epsilon::Real=1e-5) where {T<:Real, M}

Create averaged power spectrum from an event list divided into segments.

# Arguments
- `events`: EventList to analyze
- `segment_size`: Size of segments in seconds
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")
- `dt`: Bin size in seconds for creating light curves from each segment
- `epsilon`: Tolerance for segment boundaries

# Returns
- `AveragedPowerspectrum` object

# Examples
```julia
events = readevents("data.fits")
ps_avg = AveragedPowerspectrum(events, 1024.0, norm="leahy", dt=0.1)
```
"""
function AveragedPowerspectrum(events::EventList{Vector{T}, M}, segment_size::Real; 
                             norm::String="frac", 
                             dt::Real=1.0,
                             epsilon::Real=1e-5) where {T<:Real, M}
    # this segment is for debug/valid can be removed build for [Comprehensive Error Handling tests]
    if isnan(segment_size)
        throw(ArgumentError("Segment size cannot be NaN"))
    end
    if isinf(segment_size)
        throw(ArgumentError("Segment size cannot be Inf"))
    end
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive"))
    end
    # Input validation
    length(events) > 1 || throw(ArgumentError("EventList must have more than 1 event"))
    dt > 0 || throw(ArgumentError("Bin size must be positive"))
    segment_size > dt || throw(ArgumentError("Segment size must be larger than bin size"))
    
    n_bins_per_segment = round(Int, segment_size / dt)
    
    @debug "Segment parameters" segment_size dt n_bins_per_segment expected_duration=n_bins_per_segment*dt
    
    # Get GTIs from metadata
    gtis = if has_gti(events)
        events.meta.gti
    else
        # Create single GTI spanning the entire observation
        time_span = extrema(events.times)
        reshape([time_span[1], time_span[2]], 1, 2)
    end
    
    @debug "Processing EventList" total_events=length(events) segment_size=segment_size dt=dt bins_per_segment=n_bins_per_segment gti_shape=size(gtis)
    
    # unbinned segment generator for events
    segment_generator = generate_indices_of_segment_boundaries_unbinned(
        events.times, gtis, segment_size
    )
    
    # Initialize arrays
    freqs = fftfreq(n_bins_per_segment, 1/dt)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero=false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]
    
    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0
    
    # Process each segment using the generator
    for (start_time, stop_time, start_idx, stop_idx) in segment_generator
        # Extract events in this segment
        segment_events = if start_idx <= stop_idx && start_idx > 0 && stop_idx <= length(events.times)
            events.times[start_idx:stop_idx]
        else
            # Handle edge cases with time-based filtering
            filter(t -> start_time <= t < stop_time, events.times)
        end
        
        # Skip if too few events
        if length(segment_events) < 2
            @debug "Skipping segment: too few events" start_time stop_time n_events=length(segment_events)
            continue
        end
        
        # Create light curve for this segment using your existing function
        try
            # Create a temporary EventList for this segment
            segment_event_list = EventList(segment_events, nothing, events.meta)
            segment_lc = create_lightcurve(segment_event_list, dt, 
                                         tstart=start_time, tstop=stop_time)
            
            # Skip if zero counts
            segment_sum = sum(segment_lc.counts)
            if segment_sum == 0
                @debug "Skipping segment: zero counts" start_time stop_time
                continue
            end
            
            # Handle variable segment lengths by padding or truncating
            segment_counts = if length(segment_lc.counts) == n_bins_per_segment
                segment_lc.counts
            elseif length(segment_lc.counts) > n_bins_per_segment
                # Truncate if too long
                @debug "Truncating segment" actual_length=length(segment_lc.counts) expected=n_bins_per_segment
                segment_lc.counts[1:n_bins_per_segment]
            else
                # Pad with zeros if too short
                @debug "Padding segment with zeros" actual_length=length(segment_lc.counts) expected=n_bins_per_segment
                padded = zeros(eltype(segment_lc.counts), n_bins_per_segment)
                padded[1:length(segment_lc.counts)] .= segment_lc.counts
                padded
            end
            
            # Calculate FFT and power
            ft = fft(segment_lc.counts)
            unnorm_power = abs2.(ft[pos_freq_idx])
            
            # Normalize power
            power = normalize_periodograms(
                unnorm_power,
                dt,
                n_bins_per_segment;
                mean_flux = mean(segment_lc.counts),
                n_ph = segment_sum,
                norm = norm,
                power_type = "all"
            )
            
            total_power .+= power
            total_counts += segment_sum
            n_segments_used += 1
            
            @debug "Processed segment" start_time stop_time n_events=length(segment_events) n_counts=segment_sum mean_rate=segment_sum/segment_size
            
        catch e
            @debug "Failed to process segment" start_time stop_time error=e
            continue
        end
    end
    
    @debug "Segment processing complete" segments_used=n_segments_used
    
    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end
    
    # Calculate final results
    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)
    
    # Calculate errors
    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end
    
    # Create metadata for the result (your LightCurveMetadata structure)
    headers = events.meta.headers
    
    # Extract values from FITSHeader - handle both Dict-like and FITSHeader access
    telescope = try
        get(headers, "TELESCOP", "")
    catch
        haskey(headers, "TELESCOP") ? headers["TELESCOP"] : ""
    end
    
    instrument = try
        get(headers, "INSTRUME", "")
    catch
        haskey(headers, "INSTRUME") ? headers["INSTRUME"] : ""
    end
    
    object_name = try
        get(headers, "OBJECT", "")
    catch
        haskey(headers, "OBJECT") ? headers["OBJECT"] : ""
    end
    
    mjdref = try
        get(headers, "MJDREF", 0.0)
    catch
        haskey(headers, "MJDREF") ? headers["MJDREF"] : 0.0
    end
    
    result_metadata = LightCurveMetadata(
        telescope,
        instrument,
        object_name,
        mjdref,
        (minimum(events.times), maximum(events.times)),
        Float64(dt),
        [Dict{String,Any}("TELESCOP" => telescope, "INSTRUME" => instrument, "OBJECT" => object_name, "MJDREF" => mjdref)],  # Create a simple dict with key info
        Dict(
            "original_file" => events.meta.filepath,
            "original_hdu" => events.meta.hdu,
            "energy_units" => events.meta.energy_units,
            "n_original_events" => length(events),
            "gti" => gtis,
            "original_fits_header" => headers
        )
    )
    
    return AveragedPowerspectrum{T}(
        freqs,
        avg_power,
        power_err,
        norm,
        df,
        segment_size,
        total_counts,
        n_segments_used,
        mean_rate,
        length(freqs),
        result_metadata
    )
end