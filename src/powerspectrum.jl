"""
Abstract type representing a power spectrum, which characterizes the distribution 
of power across different frequencies in a signal.

Subtypes include:
- PowerSpectrum{T}: Represents a power spectrum for a single signal segment
- AveragedPowerSpectrum{T}: Represents a power spectrum averaged over multiple segments

# Type Parameters
- `T`: The numeric type for frequency and power values (typically Float64)
"""
abstract type AbstractPowerSpectrum{T} end
"""
    PowerSpectrum{T} <: AbstractPowerSpectrum{T}

Power spectrum of a single light curve segment.

# Fields
- `freq::Vector{T}`: Frequencies
- `power::Vector{T}`: Power values
- `power_err::Vector{T}`: Uncertainties on power
- `norm::String`: Normalization used
- `df::T`: Frequency resolution
- `nphots::Int`: Total number of photons
- `m::Int`: Number of segments (always 1)
- `n::Int`: Number of frequencies
- `metadata::Union{LightCurveMetadata,FITSMetadata}`: Source metadata
"""
struct PowerSpectrum{T} <: AbstractPowerSpectrum{T}
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
    metadata::Union{LightCurveMetadata,FITSMetadata}
end

"""
    AveragedPowerSpectrum{T} <: AbstractPowerSpectrum{T}

Averaged power spectrum from multiple light curve segments.

# Fields
- `freq::Vector{T}`: Frequencies
- `power::Vector{T}`: Averaged power values
- `power_err::Vector{T}`: Uncertainties on power
- `norm::String`: Normalization used
- `df::T`: Frequency resolution
- `segment_size::T`: Length of each segment in seconds
- `nphots::Int`: Total number of photons processed
- `m::Int`: Number of segments averaged
- `mean_rate::T`: Mean count rate
- `n::Int`: Number of frequencies
- `metadata::Union{LightCurveMetadata,FITSMetadata}`: Source metadata
"""
struct AveragedPowerSpectrum{T} <: AbstractPowerSpectrum{T}
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
    metadata::Union{LightCurveMetadata,FITSMetadata}
end


#lightcurve==>

"""
    PowerSpectrum(lc::LightCurve{T}; norm::String="frac") where T<:Real

Create power spectrum from a light curve.

# Arguments
- `lc`: Light curve to analyze
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")

# Returns
- `PowerSpectrum` object

# Examples
```julia
ps = PowerSpectrum(lc, norm="leahy")
```
"""
function PowerSpectrum(lc::LightCurve{T}; norm::String = "frac") where {T<:Real}
    bin_size = lc.metadata.bin_size
    n_bins = length(lc.counts)

    n_bins > 1 || throw(ArgumentError("Light curve must have more than 1 bin"))
    bin_size > 0 || throw(ArgumentError("Bin size must be positive"))

    ft = fft(lc.counts)

    # Use proper sampling frequency
    freqs = fftfreq(n_bins, 1 / bin_size)
    pos_freq_idx = positive_fft_bins(n_bins)
    freqs = freqs[pos_freq_idx]

    unnorm_power = abs2.(ft[pos_freq_idx])

    power = normalize_periodograms(
        unnorm_power,
        bin_size,
        n_bins;
        mean_flux = mean(lc.counts),
        n_ph = sum(lc.counts),
        norm = norm,
        power_type = "all",
    )

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
        lc.metadata,
    )
end

"""
    AveragedPowerSpectrum(lc::LightCurve{T}, segment_size::Real; norm::String="frac", epsilon::Real=1e-5) where T<:Real

Create averaged power spectrum from a light curve divided into segments.

# Arguments
- `lc`: Light curve to analyze
- `segment_size`: Size of segments in seconds
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")
- `epsilon`: Tolerance for segment boundaries

# Returns
- `AveragedPowerSpectrum` object

# Examples
```julia
ps_avg = AveragedPowerSpectrum(lc, 1024.0, norm="leahy")
```
"""
function AveragedPowerSpectrum(
    lc::LightCurve{T},
    segment_size::Real;
    norm::String = "frac",
    epsilon::Real = 1e-5,
) where {T<:Real}

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

    if n_bins_per_segment <= 1
        throw(ArgumentError("Segment size too small"))
    end

    # Get GTIs from metadata
    gtis = if hasfield(typeof(lc.metadata), :gti) && !isnothing(lc.metadata.gti)
        lc.metadata.gti
    elseif haskey(lc.metadata.extra, "gti")
        lc.metadata.extra["gti"]
    elseif haskey(lc.metadata.extra, "GTI")
        lc.metadata.extra["GTI"]
    else
        throw(ArgumentError("No GTI information found in metadata"))
    end

    # Use the appropriate generator based on whether data is binned
    segment_generator = generate_indices_of_segment_boundaries_binned(
        lc.time,
        gtis,
        segment_size,
        dt = bin_size,
    )

    freqs = fftfreq(n_bins_per_segment, 1 / bin_size)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero = false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]

    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0

    for (start_time, stop_time, start_idx, stop_idx) in segment_generator
        segment_length = stop_idx - start_idx
        if segment_length != n_bins_per_segment
            continue
        end

        segment_counts = @view lc.counts[start_idx+1:stop_idx]

        segment_sum = sum(segment_counts)
        if segment_sum == 0
            continue
        end

        ft = fft(segment_counts)
        unnorm_power = abs2.(ft[pos_freq_idx])

        power = normalize_periodograms(
            unnorm_power,
            bin_size,
            n_bins_per_segment;
            mean_flux = mean(segment_counts),
            n_ph = segment_sum,
            norm = norm,
            power_type = "all",
        )

        total_power .+= power
        total_counts += segment_sum
        n_segments_used += 1
    end

    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end

    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)

    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end

    return AveragedPowerSpectrum{T}(
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
        lc.metadata,
    )
end

#Eventlist==>

"""
    PowerSpectrum(events::EventList{Vector{T}, M}; norm::String="frac", dt::Real=1.0) where {T<:Real, M}

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
ps = PowerSpectrum(events, norm="leahy", dt=0.1)
```
"""
function PowerSpectrum(
    events::EventList{Vector{T},M},
    dt::Real,
    segment_size::Real;
    norm::String = "leahy",
) where {T<:Real,M}

    length(events) > 1 || throw(ArgumentError("EventList must have more than 1 event"))
    dt > 0 || throw(ArgumentError("Bin size must be positive"))
    segment_size > dt || throw(ArgumentError("Segment size must be larger than bin size"))

    n_bins_per_segment = round(Int, segment_size / dt)

    gtis = if has_gti(events)
        events.meta.gti
    else
        # Create single GTI spanning the entire observation
        time_span = extrema(events.times)
        reshape([time_span[1], time_span[2]], 1, 2)
    end

    # Use unbinned segment generator - this preserves exact event timing
    segment_generator =
        generate_indices_of_segment_boundaries_unbinned(events.times, gtis, segment_size)

    freqs = fftfreq(n_bins_per_segment, 1 / dt)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero = false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]

    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0

    for (start_time, stop_time, start_idx, stop_idx) in segment_generator

        # Extract events in this segment - preserve exact timing
        segment_event_times =
            if start_idx <= stop_idx && start_idx > 0 && stop_idx <= length(events.times)
                @view events.times[start_idx:stop_idx]
            else
                # Fallback to time-based filtering for edge cases
                filter(t -> start_time <= t < stop_time, events.times)
            end

        if length(segment_event_times) < 2
            continue
        end

        # Create time grid for this segment
        time_grid = range(start_time, stop = stop_time, length = n_bins_per_segment + 1)
        bin_centers = (time_grid[1:end-1] + time_grid[2:end]) / 2

        # Bin events directly without creating LightCurve object
        # This preserves more control over the binning process
        counts = zeros(Int, n_bins_per_segment)

        for event_time in segment_event_times
            bin_idx = searchsortedfirst(time_grid, event_time)
            if 1 <= bin_idx <= n_bins_per_segment
                counts[bin_idx] += 1
            end
        end

        segment_total_counts = sum(counts)

        if segment_total_counts == 0
            continue
        end

        ft = fft(counts)
        unnorm_power = abs2.(ft[pos_freq_idx])

        power = normalize_periodograms(
            unnorm_power,
            dt,
            n_bins_per_segment;
            mean_flux = mean(counts),
            n_ph = segment_total_counts,
            norm = norm,
            power_type = "all",
        )

        total_power .+= power
        total_counts += segment_total_counts
        n_segments_used += 1
    end

    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end

    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)

    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end

    result_metadata = create_powerspectrum_metadata(events, dt, segment_size)

    return AveragedPowerSpectrum{T}(
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
        result_metadata,
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

    gtis =
        has_gti(events) ? events.meta.gti :
        reshape([minimum(events.times), maximum(events.times)], 1, 2)

    return LightCurveMetadata(
        telescope,
        instrument,
        object_name,
        mjdref,
        (minimum(events.times), maximum(events.times)),
        Float64(dt),
        [
            Dict{String,Any}(
                "TELESCOP" => telescope,
                "INSTRUME" => instrument,
                "OBJECT" => object_name,
                "MJDREF" => mjdref,
            ),
        ],
        Dict(
            "analysis_method" => "direct_events_processing",
            "original_file" => events.meta.filepath,
            "original_hdu" => events.meta.hdu,
            "energy_units" => events.meta.energy_units,
            "n_original_events" => length(events),
            "segment_size" => segment_size,
            "time_resolution" => dt,
            "gti" => gtis,
            "original_fits_header" => headers,
        ),
    )
end
"""
    AveragedPowerSpectrum(events::EventList{Vector{T}, M}, segment_size::Real; 
                         norm::String="frac", dt::Real=1.0, 
                         epsilon::Real=1e-5) where {T<:Real, M}

Create averaged power spectrum from an event list divided into segments.
Uses direct event binning without creating intermediate LightCurve objects.

# Arguments
- `events`: EventList to analyze
- `segment_size`: Size of segments in seconds
- `norm`: Normalization type ("leahy", "frac", "rms", "abs")
- `dt`: Bin size in seconds for creating light curves from each segment
- `epsilon`: Tolerance for segment boundaries

# Returns
- `AveragedPowerSpectrum` object

# Examples
```julia
events = readevents("data.fits")
ps_avg = AveragedPowerSpectrum(events, 1024.0, norm="leahy", dt=0.1)
```
"""
function AveragedPowerSpectrum(
    events::EventList{Vector{T},M},
    segment_size::Real;
    norm::String = "frac",
    dt::Real = 1.0,
    epsilon::Real = 1e-5,
) where {T<:Real,M}

    if isnan(segment_size)
        throw(ArgumentError("Segment size cannot be NaN"))
    end
    if isinf(segment_size)
        throw(ArgumentError("Segment size cannot be Inf"))
    end
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive"))
    end
    if length(events) <= 1
        throw(ArgumentError("EventList must have more than 1 event"))
    end
    if dt <= 0
        throw(ArgumentError("Bin size must be positive"))
    end
    if segment_size <= dt
        throw(ArgumentError("Segment size must be larger than bin size"))
    end

    n_bins_per_segment = round(Int, segment_size / dt)

    if n_bins_per_segment < 2
        throw(
            ArgumentError(
                "Segment size too small relative to dt: results in < 2 bins per segment",
            ),
        )
    end

    gtis = if has_gti(events)
        events.meta.gti
    else
        time_span = extrema(events.times)
        reshape([time_span[1], time_span[2]], 1, 2)
    end

    # Use unbinned segment generator - this preserves exact event timing
    # This is the key function your mentor wants you to focus on!
    segment_generator =
        generate_indices_of_segment_boundaries_unbinned(events.times, gtis, segment_size)

    freqs = fftfreq(n_bins_per_segment, 1 / dt)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero = false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]

    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0

    for (start_time, stop_time, start_idx, stop_idx) in segment_generator

        # Extract events in this segment - preserve exact timing
        segment_event_times =
            if start_idx <= stop_idx && start_idx > 0 && stop_idx <= length(events.times)
                @view events.times[start_idx:stop_idx]
            else
                filter(t -> start_time <= t < stop_time, events.times)
            end

        if length(segment_event_times) < 2
            continue
        end

        time_grid = range(start_time, stop = stop_time, length = n_bins_per_segment + 1)

        # Bin events directly without creating LightCurve object
        # This preserves more control over the binning process
        counts = zeros(Int, n_bins_per_segment)

        for event_time in segment_event_times
            bin_idx = searchsortedfirst(time_grid, event_time)
            if 1 <= bin_idx <= n_bins_per_segment
                counts[bin_idx] += 1
            end
        end

        segment_total_counts = sum(counts)

        if segment_total_counts == 0
            continue
        end

        ft = fft(counts)
        unnorm_power = abs2.(ft[pos_freq_idx])

        power = normalize_periodograms(
            unnorm_power,
            dt,
            n_bins_per_segment;
            mean_flux = mean(counts),
            n_ph = segment_total_counts,
            norm = norm,
            power_type = "all",
        )

        total_power .+= power
        total_counts += segment_total_counts
        n_segments_used += 1
    end

    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end

    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)

    power_err = if norm == "leahy"
        fill(2.0, length(avg_power))
    elseif norm in ["frac", "rms"]
        avg_power ./ sqrt(n_segments_used)
    else
        sqrt.(avg_power ./ n_segments_used)
    end

    result_metadata = create_powerspectrum_metadata(events, dt, segment_size)

    return AveragedPowerSpectrum{T}(
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
        result_metadata,
    )
end
"""
    freqs(ps::AbstractPowerSpectrum) -> Vector{T}

Get frequency values from a power spectrum.

# Arguments
- `ps::AbstractPowerSpectrum`: Input power spectrum

# Returns
- `Vector{T}`: Frequency array in Hz

# Examples
```julia
ps = powerspectrum(events)
frequencies = freqs(ps)
```
"""
freqs(ps::AbstractPowerSpectrum) = ps.freq

"""
    power(ps::AbstractPowerSpectrum) -> Vector{T}

Get power values from a power spectrum.

# Arguments
- `ps::AbstractPowerSpectrum`: Input power spectrum

# Returns
- `Vector{T}`: Power values (units depend on normalization)

# Examples
```julia
ps = powerspectrum(events)
power_values = power(ps)
```
"""
power(ps::AbstractPowerSpectrum) = ps.power

"""
    errors(ps::AbstractPowerSpectrum) -> Vector{T}

Get error estimates for power values from a power spectrum.

# Arguments
- `ps::AbstractPowerSpectrum`: Input power spectrum

# Returns
- `Vector{T}`: Standard error estimates for power values

# Examples
```julia
ps = powerspectrum(events)
power_uncertainties = errors(ps)
```
"""
errors(ps::AbstractPowerSpectrum) = ps.power_err

# Base methods
"""
    length(ps::AbstractPowerSpectrum) -> Int

Get the number of frequency bins in the power spectrum.
"""
Base.length(ps::AbstractPowerSpectrum) = length(ps.freq)
"""
    size(ps::AbstractPowerSpectrum) -> Tuple{Int}

Get the size of the power spectrum as a tuple.
"""
Base.size(ps::AbstractPowerSpectrum) = (length(ps),)

"""
    getindex(ps::AbstractPowerSpectrum, i) -> Tuple{T, T, T}

Get the frequency, power, and error at index i.

# Returns
- `Tuple{T, T, T}`: (frequency, power, error) at index i
"""
Base.getindex(ps::AbstractPowerSpectrum, i) = (ps.freq[i], ps.power[i], ps.power_err[i])

"""
    rebin(ps::AbstractPowerSpectrum, factor::Integer)

Rebin a power spectrum by a factor `factor`.
The new frequency resolution will be `df * factor`.
The power values are averaged in the new bins.
The errors are propagated accordingly.
"""

function rebin(ps::PowerSpectrum{T}, factor::Integer) where {T}
    if factor == 1
        return ps
    end

    n_new = floor(Int, length(ps.freq) / factor)
    new_freq = zeros(T, n_new)
    new_power = zeros(T, n_new)
    new_err = zeros(T, n_new)

    for i in 1:n_new
        start_idx = (i - 1) * factor + 1
        end_idx = i * factor
        
        new_freq[i] = mean(ps.freq[start_idx:end_idx])
        new_power[i] = mean(ps.power[start_idx:end_idx])
        
        # Error propagation for the mean: sqrt(sum(sigma_i^2)) / N
        sum_sq_err = sum(abs2, ps.power_err[start_idx:end_idx])
        new_err[i] = sqrt(sum_sq_err) / factor
    end

    return PowerSpectrum{T}(
        new_freq,
        new_power,
        new_err,
        ps.norm,
        ps.df * factor,
        ps.nphots,
        ps.m,  
        n_new,
        ps.metadata
    )
end

function rebin(ps::AveragedPowerSpectrum{T}, factor::Integer) where {T}
    if factor == 1
        return ps
    end

    n_new = floor(Int, length(ps.freq) / factor)
    new_freq = zeros(T, n_new)
    new_power = zeros(T, n_new)
    new_err = zeros(T, n_new)

    for i in 1:n_new
        start_idx = (i - 1) * factor + 1
        end_idx = i * factor
        
        new_freq[i] = mean(ps.freq[start_idx:end_idx])
        new_power[i] = mean(ps.power[start_idx:end_idx])
        
        # Error propagation for the mean: sqrt(sum(sigma_i^2)) / N
        sum_sq_err = sum(abs2, ps.power_err[start_idx:end_idx])
        new_err[i] = sqrt(sum_sq_err) / factor
    end

    return AveragedPowerSpectrum{T}(
        new_freq,
        new_power,
        new_err,
        ps.norm,
        ps.df * factor,
        ps.segment_size,
        ps.nphots,
        ps.m,
        ps.mean_rate,
        n_new,
        ps.metadata
    )
end

"""
    logrebin(ps::AbstractPowerSpectrum; f::Real=0.01)

Rebin the power spectrum logarithmically.
The frequency bins increase geometrically with a factor `1 + f`.
"""
function logrebin(ps::AbstractPowerSpectrum{T}; f::Real = 0.01) where {T}
    
    min_freq = ps.freq[1]
    max_freq = ps.freq[end]
    
    # Create geometric progression of bin edges
    # start at min_freq - df/2 to capture the first bin properly
    # but for simplicity in log rebinning often we start at min_freq
    
    # Common implementation:
    # 1. generated log spaced edges
    # 2. group original bins into these new bins
    
    min_bin_width = ps.df
    
    # We need to find the edges of the new bins
    # frequency k: freq[k] = freq[0] * (1+f)^k
    
    k = 1
    current_freq = min_freq
    new_freqs = T[]
    new_powers = T[]
    new_errs = T[]
    
    # Indices in the original arrays
    start_idx = 1
    n = length(ps.freq)
    
    while start_idx <= n
        # Determine the width of the current bin in the target log grid
        # The next boundary would be current_freq * (1+f)
        # But we must ensure it's at least min_bin_width to encompass at least one original bin
        
        target_next_freq = current_freq * (1+f)
        if target_next_freq - current_freq < min_bin_width
             target_next_freq = current_freq + min_bin_width
        end
        
        # Find which original bins fall into [current_freq, target_next_freq)
        # We assume original freqs are bin centers. 
        # A simple approach: equivalent to rebinning by non-integer factors
        
        # accumulate bins until we reach the target frequency
        end_idx = start_idx
        while end_idx <= n && ps.freq[end_idx] < target_next_freq
            end_idx += 1
        end
        end_idx -= 1
        
        if end_idx < start_idx
            end_idx = start_idx
        end
        
        # Average properties in this range
        push!(new_freqs, mean(ps.freq[start_idx:end_idx]))
        push!(new_powers, mean(ps.power[start_idx:end_idx]))
        
        # Error prop assuming independent errors: sigma_mean = sqrt(sum sigma_i^2) / N
        n_averaged = end_idx - start_idx + 1
        sum_sq_err = sum(abs2, ps.power_err[start_idx:end_idx])
        push!(new_errs, sqrt(sum_sq_err) / n_averaged)
        
        current_freq = ps.freq[end_idx] # effectively moving to the next available frequency
        start_idx = end_idx + 1
    end
    
    # Construct the new object based on the type of ps
    if ps isa PowerSpectrum
        return PowerSpectrum{T}(
            new_freqs,
            new_powers,
            new_errs,
            ps.norm,
            0.0, # df is variable
            ps.nphots,
            ps.m,
            length(new_freqs),
            ps.metadata
        )
    elseif ps isa AveragedPowerSpectrum
        return AveragedPowerSpectrum{T}(
            new_freqs,
            new_powers,
            new_errs,
            ps.norm,
            0.0, # df is variable
            ps.segment_size,
            ps.nphots,
            ps.m,
            ps.mean_rate,
            length(new_freqs),
            ps.metadata
        )
    end
end
