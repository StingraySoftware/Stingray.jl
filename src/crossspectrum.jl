"""
    AbstractCrossSpectrum{T}

Abstract type for cross-spectrum analysis structures.

Type parameter `T` represents the numeric type used for calculations (typically `Float64`).
"""
abstract type AbstractCrossSpectrum{T} end

"""
    CrossSpectrum{T} <: AbstractCrossSpectrum{T}

A structure containing cross power spectrum data and associated metadata.

# Type Parameters
- `T`: Numeric type for calculations (typically `Float64`)

# Fields
$(TYPEDFIELDS)

# Examples
```julia
# Create from event lists
cs = CrossSpectrum(ev1, ev2, 100.0, 0.1)

# Create from light curves
cs = CrossSpectrum(lc1, lc2, 100.0)
```
"""
mutable struct CrossSpectrum{T} <: AbstractCrossSpectrum{T}
    "Frequencies in Hz"
    freq::Vector{T}
    "Complex cross power values"
    power::Vector{Complex{T}}
    "Power errors"
    power_err::Union{Nothing,Vector{T}}
    "Auto power spectrum of first light curve"
    ps1::Vector{T}
    "Auto power spectrum of second light curve"
    ps2::Vector{T}
    "Normalization type (leahy, frac, rms, abs, none)"
    norm::String
    "Power type (all, real, absolute)"
    power_type::String
    "Frequency resolution (Hz)"
    df::T
    "Total number of photons in first light curve"
    nphots1::T
    "Total number of photons in second light curve"
    nphots2::T
    "Number of segments (m=1: single spectrum, m>1: averaged spectrum)"
    m::Int
    "Number of frequencies"
    n::Int
    "Number of nearby averaged frequencies"
    k::Union{Int,Vector{Int}}
    "Original light curve metadata from first signal"
    metadata1::Union{LightCurveMetadata,FITSMetadata}
    "Original light curve metadata from second signal"
    metadata2::Union{LightCurveMetadata,FITSMetadata}
    "return full frequencies -ve and +ve if true else only +ve"
    fullspec::Bool
    "if set to false Independent channels else Reference band contains subject band"
    channels_overlap::Bool
    "Segment size used for averaging (nothing for single spectrum)"
    segment_size::Union{Nothing,T}
    "Mean rate of first light curve (nothing for single spectrum)"
    mean_rate1::Union{Nothing,T}
    "Mean rate of second light curve (nothing for single spectrum)"
    mean_rate2::Union{Nothing,T}
end

"""
    is_averaged(cs::CrossSpectrum) -> Bool

Check if a CrossSpectrum represents an averaged spectrum (multiple segments).

# Arguments
- `cs::CrossSpectrum`: The cross spectrum to check

# Returns
- `Bool`: `true` if averaged (m > 1), `false` if single segment

# Examples
```julia
cs_single = CrossSpectrum(lc1, lc2)
cs_averaged = AveragedCrossSpectrum(lc1, lc2, 100.0)

is_averaged(cs_single)    # false
is_averaged(cs_averaged)  # true
```
"""
is_averaged(cs::CrossSpectrum) = cs.m > 1

"""
    is_single(cs::CrossSpectrum) -> Bool

Check if a CrossSpectrum represents a single segment spectrum.

# Arguments
- `cs::CrossSpectrum`: The cross spectrum to check

# Returns
- `Bool`: `true` if single segment (m = 1), `false` if averaged

# Examples
```julia
cs_single = CrossSpectrum(lc1, lc2)
is_single(cs_single)  # true
```
"""
is_single(cs::CrossSpectrum) = cs.m == 1

"""
    CrossSpectrum(ev1::EventList, ev2::EventList, segment_size::Real, dt::Real;
                  norm::String="frac", use_common_mean::Bool=true,
                  fullspec::Bool=false, power_type::String="all") -> CrossSpectrum

Create a single-segment cross spectrum from two event lists.

# Arguments
- `ev1::EventList`: First event list
- `ev2::EventList`: Second event list
- `segment_size::Real`: Size of the segment in time units
- `dt::Real`: Time bin size

# Keywords
- `norm::String="frac"`: Normalization type ("frac", "leahy", "rms", "abs", "none")
- `use_common_mean::Bool=true`: Whether to use common mean for normalization
- `fullspec::Bool=false`: If true, include negative frequencies
- `power_type::String="all"`: Power type ("all", "real", "absolute")

# Returns
- `CrossSpectrum`: Single-segment cross spectrum

# Examples
```julia
# Basic usage
cs = CrossSpectrum(ev1, ev2, 100.0, 0.1)

# With different normalization
cs = CrossSpectrum(ev1, ev2, 100.0, 0.1, norm="leahy")
```
"""
function CrossSpectrum(
    ev1::EventList,
    ev2::EventList,
    segment_size::Real,
    dt::Real;
    norm::String = "frac",
    use_common_mean::Bool = true,
    fullspec::Bool = false,
    power_type::String = "all",
)

    validate_time_alignment(ev1.meta.headers, ev2.meta.headers)

    validate_time_alignment(ev1.meta.headers, ev2.meta.headers)

    # Get GTIs and find intersection
    gti1 = ev1.meta.gti
    gti2 = ev2.meta.gti
    cross_gti = intersect_gtis(gti1, gti2)

    if size(cross_gti, 1) == 0
        error("No overlapping GTIs between event lists!")
    end

    # Calculate number of bins
    n_bin = round(Int, segment_size / dt)
    adjusted_dt = segment_size / n_bin

    # Create frequency array
    freq = fftfreq(n_bin, 1 / adjusted_dt)
    if !fullspec
        fgt0 = positive_fft_bins(n_bin)
        freq = freq[fgt0]
    end

    # Find and process ONLY THE FIRST valid segment
    cross_power = nothing
    pds1_power = nothing
    pds2_power = nothing
    segment_photons1 = 0.0
    segment_photons2 = 0.0

    for (s, e, idx0_1, idx1_1) in
        generate_indices_of_segment_boundaries_unbinned(ev1.times, cross_gti, segment_size)
        # Find corresponding indices in second event list
        idx0_2 = searchsortedfirst(ev2.times, s)
        idx1_2 = searchsortedfirst(ev2.times, e)

        if idx1_1 - idx0_1 < 2 || idx1_2 - idx0_2 < 2
            continue
        end

        # Get event times for this segment
        segment_times1 = @view ev1.times[idx0_1:idx1_1-1]
        segment_times2 = @view ev2.times[idx0_2:idx1_2-1]

        # Create binned light curves
        edges = range(s, stop = e, length = n_bin + 1)
        lc1 = fit(Histogram, segment_times1, edges).weights
        lc2 = fit(Histogram, segment_times2, edges).weights

        if length(lc1) != n_bin || length(lc2) != n_bin
            continue
        end

        # Calculate FFTs
        ft1 = fft(Float64.(lc1))
        ft2 = fft(Float64.(lc2))

        # Calculate cross spectrum and power spectra
        cross_power = ft1 .* conj.(ft2)
        pds1_power = abs2.(ft1)
        pds2_power = abs2.(ft2)

        # Take only positive frequencies if needed
        if !fullspec
            cross_power = cross_power[fgt0]
            pds1_power = pds1_power[fgt0]
            pds2_power = pds2_power[fgt0]
        end

        # Store photon counts
        segment_photons1 = sum(lc1)
        segment_photons2 = sum(lc2)

        # BREAK after first valid segment
        break
    end

    if isnothing(cross_power)
        error("No valid segments found for cross spectrum")
    end

    # Calculate mean rates for this single segment
    mean_rate1 = segment_photons1 / segment_size
    mean_rate2 = segment_photons2 / segment_size

    # Normalize the spectra
    if norm != "none"
        cross_power = normalize_periodograms(
            cross_power,
            adjusted_dt,
            n_bin;
            mean_flux = sqrt(segment_photons1 * segment_photons2) / n_bin,
            n_ph = sqrt(segment_photons1 * segment_photons2),
            norm = norm,
            power_type = power_type,
        )
        pds1_power = normalize_periodograms(
            pds1_power,
            adjusted_dt,
            n_bin;
            mean_flux = segment_photons1 / n_bin,
            n_ph = segment_photons1,
            norm = norm,
            power_type = power_type,
        )
        pds2_power = normalize_periodograms(
            pds2_power,
            adjusted_dt,
            n_bin;
            mean_flux = segment_photons2 / n_bin,
            n_ph = segment_photons2,
            norm = norm,
            power_type = power_type,
        )
    end

    return CrossSpectrum{Float64}(
        freq,
        cross_power,
        nothing,  # power_err
        pds1_power,
        pds2_power,
        norm,
        power_type,
        freq[2] - freq[1],  # df
        segment_photons1,
        segment_photons2,
        1,  # m = 1 for SINGLE segment
        n_bin,       # n
        1,           # k
        ev1.meta,    # metadata1
        ev2.meta,    # metadata2
        fullspec,
        false,       # channels_overlap
        Float64(segment_size),  # segment_size
        mean_rate1,  # mean_rate1
        mean_rate2,   # mean_rate2
    )
end

"""
    CrossSpectrum(lc1::LightCurve, lc2::LightCurve, segment_size::Union{Nothing,Real}=nothing;
                  norm::String="frac", use_common_mean::Bool=true,
                  fullspec::Bool=false, power_type::String="all") -> CrossSpectrum

Create a single-segment cross spectrum from two light curves.

# Arguments
- `lc1::LightCurve`: First light curve
- `lc2::LightCurve`: Second light curve  
- `segment_size::Union{Nothing,Real}=nothing`: Segment size (if nothing, uses minimum segment size)

# Keywords
- `norm::String="frac"`: Normalization type
- `use_common_mean::Bool=true`: Whether to use common mean for normalization
- `fullspec::Bool=false`: If true, include negative frequencies
- `power_type::String="all"`: Power type

# Returns
- `CrossSpectrum`: Single-segment cross spectrum

# Examples
```julia
# Auto-determine segment size
cs = CrossSpectrum(lc1, lc2)

# Specify segment size
cs = CrossSpectrum(lc1, lc2, 100.0)
```
"""
function CrossSpectrum(
    lc1::LightCurve,
    lc2::LightCurve,
    segment_size::Union{Nothing,Real} = nothing;
    norm::String = "frac",
    use_common_mean::Bool = true,
    fullspec::Bool = false,
    power_type::String = "all",
)

    gti1 = haskey(lc1.metadata.extra, "gti") ? lc1.metadata.extra["gti"] : nothing
    gti2 = haskey(lc2.metadata.extra, "gti") ? lc2.metadata.extra["gti"] : nothing
    cross_gti =
        isnothing(gti1) ? gti2 : (isnothing(gti2) ? gti1 : intersect_gtis(gti1, gti2))

    if isnothing(cross_gti) || size(cross_gti, 1) == 0
        error("No overlapping GTIs between light curves!")
    end

    # Apply GTIs to get segments
    lcs1 = apply_gtis(lc1, cross_gti)
    lcs2 = apply_gtis(lc2, cross_gti)
    minseg = min(length(lcs1), length(lcs2))

    if minseg == 0
        error("No valid segments after GTI intersection!")
    end

    # Get time bin size
    dt_val = isa(lc1.dt, Vector) ? lc1.dt[1] : lc1.dt

    if isnothing(segment_size)
        segment_size = minimum([lcs1[i].dt * length(lcs1[i].counts) for i = 1:minseg])
    end

    # Calculate number of bins
    n_bin = round(Int, segment_size / dt_val)
    adjusted_dt = segment_size / n_bin

    # Create frequency array
    freq = fftfreq(n_bin, 1 / adjusted_dt)
    if !fullspec
        fgt0 = positive_fft_bins(n_bin)
        freq = freq[fgt0]
    end

    # Process ONLY THE FIRST valid segment
    cross_power = nothing
    pds1_power = nothing
    pds2_power = nothing
    segment_photons1 = 0.0
    segment_photons2 = 0.0

    for i = 1:minseg
        lc1_seg = lcs1[i]
        lc2_seg = lcs2[i]
        
        # Handle different segment lengths by truncating to common length
        min_len = min(length(lc1_seg.counts), length(lc2_seg.counts))
        
        if min_len < 2
            continue
        end
        
        # Use common length
        counts1 = lc1_seg.counts[1:min_len]
        counts2 = lc2_seg.counts[1:min_len]
        
        # Ensure we have the expected number of bins or adjust
        if length(counts1) != n_bin || length(counts2) != n_bin
            # If segment is longer than expected, truncate
            if length(counts1) >= n_bin && length(counts2) >= n_bin
                counts1 = counts1[1:n_bin]
                counts2 = counts2[1:n_bin]
            else
                # If too short, skip this segment
                continue
            end
        end

        # Calculate FFTs
        ft1 = fft(Float64.(counts1))
        ft2 = fft(Float64.(counts2))

        # Calculate cross spectrum and power spectra
        cross_power = ft1 .* conj.(ft2)
        pds1_power = abs2.(ft1)
        pds2_power = abs2.(ft2)

        # Take only positive frequencies if needed
        if !fullspec
            cross_power = cross_power[fgt0]
            pds1_power = pds1_power[fgt0]
            pds2_power = pds2_power[fgt0]
        end

        # Store photon counts
        segment_photons1 = sum(counts1)
        segment_photons2 = sum(counts2)

        # BREAK after first valid segment - this is the key difference!
        break
    end

    if isnothing(cross_power)
        error("No valid segments found for cross spectrum")
    end

    # Calculate mean rates for this single segment
    mean_rate1 = segment_photons1 / segment_size
    mean_rate2 = segment_photons2 / segment_size

    # Normalize the spectra
    if norm != "none"
        cross_power = normalize_periodograms(
            cross_power,
            adjusted_dt,
            n_bin;
            mean_flux = sqrt(segment_photons1 * segment_photons2) / n_bin,
            n_ph = sqrt(segment_photons1 * segment_photons2),
            norm = norm,
            power_type = power_type,
        )
        pds1_power = normalize_periodograms(
            pds1_power,
            adjusted_dt,
            n_bin;
            mean_flux = segment_photons1 / n_bin,
            n_ph = segment_photons1,
            norm = norm,
            power_type = power_type,
        )
        pds2_power = normalize_periodograms(
            pds2_power,
            adjusted_dt,
            n_bin;
            mean_flux = segment_photons2 / n_bin,
            n_ph = segment_photons2,
            norm = norm,
            power_type = power_type,
        )
    end

    return CrossSpectrum{Float64}(
        freq,
        cross_power,
        nothing,  # power_err
        pds1_power,
        pds2_power,
        norm,
        power_type,
        freq[2] - freq[1],  # df
        segment_photons1,
        segment_photons2,
        1,  # m = 1 for SINGLE segment
        n_bin,       # n
        1,           # k
        lc1.metadata,  # metadata1
        lc2.metadata,  # metadata2
        fullspec,
        false,          # channels_overlap
        Float64(segment_size),  # segment_size
        mean_rate1,     # mean_rate1
        mean_rate2,      # mean_rate2
    )
end

"""
    AveragedCrossSpectrum(lc1::LightCurve{T}, lc2::LightCurve{T}, segment_size::Real; 
                          norm::String="frac", use_common_mean::Bool=true,
                          fullspec::Bool=false, power_type::String="all",
                          fill_errors_on_creation::Bool=true) where T<:Real -> CrossSpectrum

Create an averaged cross spectrum from two light curves by averaging multiple segments.

# Arguments
- `lc1::LightCurve{T}`: First light curve
- `lc2::LightCurve{T}`: Second light curve
- `segment_size::Real`: Size of each segment for averaging

# Keywords
- `norm::String="frac"`: Normalization type
- `use_common_mean::Bool=true`: Whether to use common mean for normalization
- `fullspec::Bool=false`: If true, include negative frequencies
- `power_type::String="all"`: Power type
- `fill_errors_on_creation::Bool=true`: Whether to calculate errors immediately

# Returns
- `CrossSpectrum`: Multi-segment averaged cross spectrum with m > 1

# Examples
```julia
# Create averaged cross spectrum
cs_avg = AveragedCrossSpectrum(lc1, lc2, 100.0)

# Without error calculation
cs_avg = AveragedCrossSpectrum(lc1, lc2, 100.0, fill_errors_on_creation=false)
```
"""
function AveragedCrossSpectrum(
    lc1::LightCurve{T},
    lc2::LightCurve{T},
    segment_size::Real;
    norm::String = "frac",
    use_common_mean::Bool = true,
    fullspec::Bool = false,
    power_type::String = "all",
    fill_errors_on_creation::Bool = true,
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

    bin_size1 = lc1.metadata.bin_size
    bin_size2 = lc2.metadata.bin_size

    if abs(bin_size1 - bin_size2) > 1e-10
        throw(ArgumentError("Light curves must have the same bin size"))
    end

    bin_size = bin_size1
    n_bins_per_segment = round(Int, segment_size / bin_size)

    if n_bins_per_segment <= 1
        throw(ArgumentError("Segment size too small"))
    end

    # Extract GTI information
    gti1 = if hasfield(typeof(lc1.metadata), :gti) && !isnothing(lc1.metadata.gti)
        lc1.metadata.gti
    elseif haskey(lc1.metadata.extra, "gti")
        lc1.metadata.extra["gti"]
    elseif haskey(lc1.metadata.extra, "GTI")
        lc1.metadata.extra["GTI"]
    else
        throw(ArgumentError("No GTI information found in first light curve metadata"))
    end

    gti2 = if hasfield(typeof(lc2.metadata), :gti) && !isnothing(lc2.metadata.gti)
        lc2.metadata.gti
    elseif haskey(lc2.metadata.extra, "gti")
        lc2.metadata.extra["gti"]
    elseif haskey(lc2.metadata.extra, "GTI")
        lc2.metadata.extra["GTI"]
    else
        throw(ArgumentError("No GTI information found in second light curve metadata"))
    end

    cross_gti = intersect_gtis(gti1, gti2)

    if size(cross_gti, 1) == 0
        throw(ArgumentError("No overlapping GTIs between light curves"))
    end

    # Generate segment boundaries
    segment_generator1 = generate_indices_of_segment_boundaries_binned(
        lc1.time,
        cross_gti,
        segment_size,
        dt = bin_size,
    )
    segment_generator2 = generate_indices_of_segment_boundaries_binned(
        lc2.time,
        cross_gti,
        segment_size,
        dt = bin_size,
    )

    # Create frequency array
    freqs = fftfreq(n_bins_per_segment, 1 / bin_size)
    if !fullspec
        pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero = false)
        freqs = freqs[pos_freq_idx]
    else
        pos_freq_idx = 1:length(freqs)
    end
    df = freqs[2] - freqs[1]

    # Initialize accumulators
    total_cross_power = zeros(Complex{T}, length(pos_freq_idx))
    total_power1 = zeros(T, length(pos_freq_idx))
    total_power2 = zeros(T, length(pos_freq_idx))
    total_counts1 = 0
    total_counts2 = 0
    n_segments_used = 0

    segments1 = collect(segment_generator1)
    segments2 = collect(segment_generator2)

    if length(segments1) != length(segments2)
        throw(ArgumentError("Mismatch in number of segments between light curves"))
    end

    # Process each segment
    for i = 1:length(segments1)
        start_time1, stop_time1, start_idx1, stop_idx1 = segments1[i]
        start_time2, stop_time2, start_idx2, stop_idx2 = segments2[i]

        segment_length1 = stop_idx1 - start_idx1
        segment_length2 = stop_idx2 - start_idx2

        if segment_length1 != n_bins_per_segment || segment_length2 != n_bins_per_segment
            continue
        end

        segment_counts1 = @view lc1.counts[start_idx1+1:stop_idx1]
        segment_counts2 = @view lc2.counts[start_idx2+1:stop_idx2]

        segment_sum1 = sum(segment_counts1)
        segment_sum2 = sum(segment_counts2)

        if segment_sum1 == 0 || segment_sum2 == 0
            continue
        end

        # Calculate FFTs
        ft1 = fft(segment_counts1)
        ft2 = fft(segment_counts2)

        unnorm_cross_power = ft1 .* conj.(ft2)
        unnorm_power1 = abs2.(ft1)
        unnorm_power2 = abs2.(ft2)

        if !fullspec
            unnorm_cross_power = unnorm_cross_power[pos_freq_idx]
            unnorm_power1 = unnorm_power1[pos_freq_idx]
            unnorm_power2 = unnorm_power2[pos_freq_idx]
        end

        # Normalize
        cross_power = normalize_periodograms(
            unnorm_cross_power,
            bin_size,
            n_bins_per_segment;
            mean_flux = sqrt(mean(segment_counts1) * mean(segment_counts2)),
            n_ph = sqrt(segment_sum1 * segment_sum2),
            norm = norm,
            power_type = power_type,
        )

        power1 = normalize_periodograms(
            unnorm_power1,
            bin_size,
            n_bins_per_segment;
            mean_flux = mean(segment_counts1),
            n_ph = segment_sum1,
            norm = norm,
            power_type = power_type,
        )

        power2 = normalize_periodograms(
            unnorm_power2,
            bin_size,
            n_bins_per_segment;
            mean_flux = mean(segment_counts2),
            n_ph = segment_sum2,
            norm = norm,
            power_type = power_type,
        )

        # Accumulate
        total_cross_power .+= cross_power
        total_power1 .+= power1
        total_power2 .+= power2
        total_counts1 += segment_sum1
        total_counts2 += segment_sum2
        n_segments_used += 1
    end

    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end

    # Average
    avg_cross_power = total_cross_power ./ n_segments_used
    avg_power1 = total_power1 ./ n_segments_used
    avg_power2 = total_power2 ./ n_segments_used

    # Calculate mean rates
    mean_rate1 = total_counts1 / (n_segments_used * segment_size)
    mean_rate2 = total_counts2 / (n_segments_used * segment_size)

    # Create object using unified CrossSpectrum struct
    cs = CrossSpectrum{T}(
        freqs,
        avg_cross_power,
        nothing,  # Will be filled by fill_errors! if requested
        avg_power1,
        avg_power2,
        norm,
        power_type,
        df,
        T(total_counts1),
        T(total_counts2),
        n_segments_used,  # m > 1 indicates averaged
        length(freqs),
        1,
        lc1.metadata,
        lc2.metadata,
        fullspec,
        false,
        T(segment_size),
        mean_rate1,
        mean_rate2,
    )

    # Fill proper errors if requested
    if fill_errors_on_creation
        fill_errors!(cs)
    end

    return cs
end

"""
    AveragedCrossSpectrum(ev1::EventList, ev2::EventList, segment_size::Real, dt::Real;
                          norm::String="frac", use_common_mean::Bool=true,
                          fullspec::Bool=false, power_type::String="all",
                          fill_errors_on_creation::Bool=true) -> CrossSpectrum

Create an averaged cross spectrum from two event lists by averaging multiple segments.

# Arguments
- `ev1::EventList`: First event list
- `ev2::EventList`: Second event list
- `segment_size::Real`: Size of each segment for averaging
- `dt::Real`: Time bin size

# Keywords
- `norm::String="frac"`: Normalization type
- `use_common_mean::Bool=true`: Whether to use common mean for normalization
- `fullspec::Bool=false`: If true, include negative frequencies
- `power_type::String="all"`: Power type
- `fill_errors_on_creation::Bool=true`: Whether to calculate errors immediately

# Returns
- `CrossSpectrum`: Multi-segment averaged cross spectrum

# Examples
```julia
# Create averaged cross spectrum from events
cs_avg = AveragedCrossSpectrum(ev1, ev2, 100.0, 0.1)
```
"""

function AveragedCrossSpectrum(
    ev1::EventList,
    ev2::EventList,
    segment_size::Real,
    dt::Real;
    norm::String = "frac",
    use_common_mean::Bool = true,
    fullspec::Bool = false,
    power_type::String = "all",
    fill_errors_on_creation::Bool = true,
)

    if isnan(segment_size)
        throw(ArgumentError("Segment size cannot be NaN"))
    end
    if isinf(segment_size)
        throw(ArgumentError("Segment size cannot be Inf"))
    end
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive"))
    end
    if dt <= 0
        throw(ArgumentError("Bin size must be positive"))
    end
    if segment_size <= dt
        throw(ArgumentError("Segment size must be larger than bin size"))
    end

    # Validate time alignment
    validate_time_alignment(ev1.meta.headers, ev2.meta.headers)

    # Get GTI intersection
    gti1 = ev1.meta.gti
    gti2 = ev2.meta.gti
    cross_gti = intersect_gtis(gti1, gti2)

    if size(cross_gti, 1) == 0
        throw(ArgumentError("No overlapping GTIs between event lists"))
    end

    n_bins_per_segment = round(Int, segment_size / dt)

    if n_bins_per_segment < 2
        throw(ArgumentError("Segment size too small relative to dt"))
    end

    # Generate segments
    segment_generator =
        generate_indices_of_segment_boundaries_unbinned(ev1.times, cross_gti, segment_size)

    # Create frequency array
    freqs = fftfreq(n_bins_per_segment, 1 / dt)
    if !fullspec
        pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero = false)
        freqs = freqs[pos_freq_idx]
    else
        pos_freq_idx = 1:length(freqs)
    end
    df = freqs[2] - freqs[1]

    # Initialize accumulators
    total_cross_power = zeros(Complex{Float64}, length(pos_freq_idx))
    total_power1 = zeros(Float64, length(pos_freq_idx))
    total_power2 = zeros(Float64, length(pos_freq_idx))
    total_counts1 = 0.0
    total_counts2 = 0.0
    n_segments_used = 0

    # Process segments
    for (start_time, stop_time, start_idx1, stop_idx1) in segment_generator
        start_idx2 = searchsortedfirst(ev2.times, start_time)
        stop_idx2 = searchsortedfirst(ev2.times, stop_time)

        if stop_idx1 - start_idx1 < 2 || stop_idx2 - start_idx2 < 2
            continue
        end

        segment_times1 = @view ev1.times[start_idx1:stop_idx1-1]
        segment_times2 = @view ev2.times[start_idx2:stop_idx2-1]

        # Create time grid and bin events
        time_grid = range(start_time, stop = stop_time, length = n_bins_per_segment + 1)

        counts1 = zeros(Int, n_bins_per_segment)
        counts2 = zeros(Int, n_bins_per_segment)

        for event_time in segment_times1
            bin_idx = searchsortedfirst(time_grid, event_time)
            if 1 <= bin_idx <= n_bins_per_segment
                counts1[bin_idx] += 1
            end
        end

        for event_time in segment_times2
            bin_idx = searchsortedfirst(time_grid, event_time)
            if 1 <= bin_idx <= n_bins_per_segment
                counts2[bin_idx] += 1
            end
        end

        segment_sum1 = sum(counts1)
        segment_sum2 = sum(counts2)

        if segment_sum1 == 0 || segment_sum2 == 0
            continue
        end

        # Calculate FFTs
        ft1 = fft(counts1)
        ft2 = fft(counts2)

        unnorm_cross_power = ft1 .* conj.(ft2)
        unnorm_power1 = abs2.(ft1)
        unnorm_power2 = abs2.(ft2)

        if !fullspec
            unnorm_cross_power = unnorm_cross_power[pos_freq_idx]
            unnorm_power1 = unnorm_power1[pos_freq_idx]
            unnorm_power2 = unnorm_power2[pos_freq_idx]
        end

        # Normalize
        cross_power = normalize_periodograms(
            unnorm_cross_power,
            dt,
            n_bins_per_segment;
            mean_flux = sqrt(mean(counts1) * mean(counts2)),
            n_ph = sqrt(segment_sum1 * segment_sum2),
            norm = norm,
            power_type = power_type,
        )

        power1 = normalize_periodograms(
            unnorm_power1,
            dt,
            n_bins_per_segment;
            mean_flux = mean(counts1),
            n_ph = segment_sum1,
            norm = norm,
            power_type = power_type,
        )

        power2 = normalize_periodograms(
            unnorm_power2,
            dt,
            n_bins_per_segment;
            mean_flux = mean(counts2),
            n_ph = segment_sum2,
            norm = norm,
            power_type = power_type,
        )

        # Accumulate
        total_cross_power .+= cross_power
        total_power1 .+= power1
        total_power2 .+= power2
        total_counts1 += segment_sum1
        total_counts2 += segment_sum2
        n_segments_used += 1
    end

    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end

    # Average
    avg_cross_power = total_cross_power ./ n_segments_used
    avg_power1 = total_power1 ./ n_segments_used
    avg_power2 = total_power2 ./ n_segments_used

    # Calculate mean rates
    mean_rate1 = total_counts1 / (n_segments_used * segment_size)
    mean_rate2 = total_counts2 / (n_segments_used * segment_size)

    # Create object using unified CrossSpectrum struct
    cs = CrossSpectrum{Float64}(
        freqs,
        avg_cross_power,
        nothing,  # Will be filled by fill_errors! if requested
        avg_power1,
        avg_power2,
        norm,
        power_type,
        df,
        total_counts1,
        total_counts2,
        n_segments_used,  # m > 1 indicates averaged
        length(freqs),
        1,
        ev1.meta,
        ev2.meta,
        fullspec,
        false,
        Float64(segment_size),
        mean_rate1,
        mean_rate2,
    )

    # Fill proper errors if requested
    if fill_errors_on_creation
        fill_errors!(cs)
    end

    return cs
end
"""
    theoretical_noise_level(cs::CrossSpectrum{T}) where T -> T

Calculate theoretical noise level for cross spectrum based on Poisson statistics.

For cross spectra, the noise level is calculated as the geometric mean of individual 
power spectral noise levels. For averaged spectra, the noise is reduced by √m where 
m is the number of segments averaged.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `T`: Theoretical cross spectrum noise level

# Mathematical Details
- For single spectra: σ_cross = √(noise1 × noise2)
- For averaged spectra: σ_cross = √(noise1 × noise2) / √m
- Individual noise levels computed using Poisson statistics

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
noise = theoretical_noise_level(cs)
println("Theoretical noise level: ", noise)
```
"""
function theoretical_noise_level(cs::CrossSpectrum{T}) where {T}

    # For cross spectra, noise is related to individual power spectra
    # The cross spectrum noise level should be approximately:
    # σ_cross ≈ sqrt(P1 * P2) / sqrt(m) for signal-dominated regions
    # σ_cross ≈ sqrt(noise1 * noise2) / sqrt(m) for noise-dominated regions

    # Individual Poisson noise levels (per segment)
    if is_averaged(cs)
        noise1 = poisson_level(cs.norm; meanrate = cs.mean_rate1, n_ph = cs.nphots1 / cs.m)
        noise2 = poisson_level(cs.norm; meanrate = cs.mean_rate2, n_ph = cs.nphots2 / cs.m)
        # For cross spectrum, theoretical noise is the geometric mean of individual noise levels
        # divided by sqrt of number of segments averaged
        cross_noise = sqrt(noise1 * noise2) / sqrt(cs.m)
    else
        noise1 = poisson_level(cs.norm; meanrate = cs.mean_rate1, n_ph = cs.nphots1)
        noise2 = poisson_level(cs.norm; meanrate = cs.mean_rate2, n_ph = cs.nphots2)
        cross_noise = sqrt(noise1 * noise2)
    end

    return cross_noise
end

"""
    fill_errors!(cs::CrossSpectrum{T}) where T -> CrossSpectrum{T}

Fill proper error estimates for cross spectrum power values.

Computes error estimates based on theoretical noise levels and, for averaged spectra,
includes sample variance contributions where signal-to-noise ratio exceeds 1.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to fill errors for (modified in-place)

# Returns
- `CrossSpectrum{T}`: The same cross spectrum with `power_err` field populated

# Mathematical Details
- Base errors start with theoretical noise level
- For averaged spectra with S/N > 1: σ_total = √(σ_noise² + σ_sample²)
- Sample variance scales as signal/√m for averaged spectra

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
fill_errors!(cs)
# cs.power_err now contains error estimates
```

# Notes
This function modifies the input cross spectrum in-place by setting the `power_err` field.
Work is still in progress for optimal error estimation methods.
"""
function fill_errors!(cs::CrossSpectrum{T}) where {T}

    # Get theoretical noise level
    noise_level = theoretical_noise_level(cs)

    # Initialize with theoretical noise level
    errors = fill(noise_level, length(cs.power))

    if is_averaged(cs)
        # For averaged cross spectra, add sample variance contribution
        signal_power = abs.(cs.power)
        snr = signal_power ./ noise_level

        # Where S/N > 1, add sample variance term
        high_snr_mask = snr .> 1.0
        if any(high_snr_mask)
            # Sample variance scales as signal/sqrt(m)
            sample_variance = signal_power[high_snr_mask] ./ sqrt(cs.m)
            errors[high_snr_mask] =
                sqrt.(errors[high_snr_mask] .^ 2 .+ sample_variance .^ 2)
        end
    else
        # For single spectra, errors are primarily Poisson noise
        # but we can add additional terms if needed
        # errors already initialized with theoretical noise level
    end

    cs.power_err = errors
    return cs
end

"""
    white_noise_level(cs::CrossSpectrum{T}; high_freq_fraction::Real=0.2) where T -> T

Estimate white noise level from high frequency portion of spectrum.

Uses the highest frequency bins to empirically estimate the noise floor, assuming
that high frequencies are dominated by white noise rather than signal.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze
- `high_freq_fraction::Real=0.2`: Fraction of highest frequencies to use for estimation

# Returns
- `T`: Estimated white noise level

# Mathematical Details
- Uses median of high frequency power values to avoid outliers
- Compares with theoretical noise level as sanity check
- Warns if estimated noise significantly exceeds theoretical expectations

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
empirical_noise = white_noise_level(cs, high_freq_fraction=0.3)
theoretical_noise = theoretical_noise_level(cs)
println("Empirical vs theoretical: ", empirical_noise, " vs ", theoretical_noise)
```

# Warnings
- Issues warning if fewer than 5 high frequency bins available
- Issues warning if estimated noise > 3× theoretical noise
"""
function white_noise_level(cs::CrossSpectrum{T}; high_freq_fraction::Real = 0.2) where {T}

    # Use the highest frequencies to estimate noise floor
    freq_cutoff = maximum(cs.freq) * (1.0 - high_freq_fraction)
    mask = cs.freq .>= freq_cutoff

    if sum(mask) < 5
        @warn "Too few high frequency points for reliable noise estimation"
        return theoretical_noise_level(cs)
    end

    high_freq_power = abs.(cs.power[mask])

    # Use median to avoid outliers
    estimated_noise = median(high_freq_power)

    # Compare with theoretical expectation
    theoretical_noise = theoretical_noise_level(cs)

    # If estimated noise is much larger than theoretical, might indicate problems
    if estimated_noise > 3.0 * theoretical_noise
        @warn "Estimated noise level ($(estimated_noise)) much higher than theoretical ($(theoretical_noise))"
    end

    return estimated_noise
end

"""
    noise_corrected_power(cs::CrossSpectrum{T}) where T -> Vector{T}

Return power spectrum with noise floor subtracted.

Subtracts the empirically estimated white noise level from the power spectrum
amplitude, ensuring non-negative results.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to correct

# Returns
- `Vector{T}`: Noise-corrected power amplitudes

# Mathematical Details
- Uses `white_noise_level()` to estimate noise floor
- Corrected power = max(|power| - noise_level, 0.01 × noise_level)
- Minimum floor prevents negative or zero values

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
raw_power = abs.(cs.power)
corrected_power = noise_corrected_power(cs)
plot(cs.freq, [raw_power corrected_power], labels=["Raw" "Corrected"])
```
"""
function noise_corrected_power(cs::CrossSpectrum{T}) where {T}

    # Use white noise level from high frequencies
    noise_level = white_noise_level(cs)
    power_amplitude = abs.(cs.power)

    # Subtract noise floor, ensuring non-negative result
    corrected_power = max.(power_amplitude .- noise_level, 0.01 * noise_level)

    return corrected_power
end

"""
    signal_to_noise_ratio(cs::CrossSpectrum{T}) where T -> Vector{T}

Calculate signal-to-noise ratio for each frequency bin.

Computes the ratio of signal power to noise level using empirical noise
estimation from high frequencies.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `Vector{T}`: Signal-to-noise ratio for each frequency bin

# Mathematical Details
- S/N = |power(f)| / noise_level
- Uses `white_noise_level()` for noise estimation
- Higher values indicate stronger signal relative to noise

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
snr = signal_to_noise_ratio(cs)
significant_mask = snr .> 3.0
println("Significant detections: ", sum(significant_mask))
```
"""
function signal_to_noise_ratio(cs::CrossSpectrum{T}) where {T}
    # Use empirical noise estimate from high frequencies
    noise_level = white_noise_level(cs)
    signal_power = abs.(cs.power)

    return signal_power ./ noise_level
end

"""
    detect_aliasing(cs::CrossSpectrum{T}) where T -> (Bool, String)

Detect potential aliasing in cross spectrum.

Analyzes the spectrum for signs of aliasing by comparing power levels
at high and low frequencies. Aliasing typically manifests as artificially
high power at frequencies near the Nyquist frequency.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `Tuple{Bool, String}`: (aliasing_detected, diagnostic_message)

# Mathematical Details
- Compares median power in highest 20% vs lowest 20% of frequencies
- Aliasing suspected if: high_freq_power > 5× low_freq_power AND > 10× noise_level
- Uses theoretical noise level for comparison

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
aliased, message = detect_aliasing(cs)
if aliased
    println("Warning: ", message)
    # Consider using anti-aliasing filters or higher sampling rate
end
```

# Detection Criteria
- High frequency power significantly exceeds low frequency power
- High frequency power significantly exceeds expected noise level
- Requires at least 5 bins in both high and low frequency regions
"""
function detect_aliasing(cs::CrossSpectrum{T}) where {T}

    # Check if power increases dramatically at high frequencies
    # This is often a sign of aliasing

    nyquist_freq = maximum(cs.freq)
    high_freq_mask = cs.freq .> 0.8 * nyquist_freq
    low_freq_mask = cs.freq .< 0.2 * nyquist_freq

    if sum(high_freq_mask) < 5 || sum(low_freq_mask) < 5
        return false, "Not enough frequency bins for aliasing detection"
    end

    high_freq_power = median(abs.(cs.power[high_freq_mask]))
    low_freq_power = median(abs.(cs.power[low_freq_mask]))

    # If high frequency power is much larger than low frequency power,
    # and much larger than expected noise, suspect aliasing
    noise_level = theoretical_noise_level(cs)

    aliasing_suspected =
        (high_freq_power > 5.0 * low_freq_power) && (high_freq_power > 10.0 * noise_level)

    message = if aliasing_suspected
        "Possible aliasing detected: high freq power = $(high_freq_power), low freq power = $(low_freq_power)"
    else
        "No obvious aliasing detected"
    end

    return aliasing_suspected, message
end

"""
    coherence(cs::CrossSpectrum{T}) where T -> Vector{T}

Calculate coherence function for cross spectrum.

Coherence measures the linear correlation between two signals as a function
of frequency, ranging from 0 (no correlation) to 1 (perfect correlation).

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `Vector{T}`: Coherence values for each frequency bin

# Mathematical Details
- For averaged spectra: γ²(f) = |P_xy(f)|² / (P_xx(f) × P_yy(f))
- For single spectra: Coherence is theoretically 1.0 everywhere
- Bias correction applied: subtracts 1/m for averaged spectra
- Values clamped to [0, 1] range

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
coh = coherence(cs)
plot(cs.freq, coh, xlabel="Frequency", ylabel="Coherence")
```

# Notes
- Only meaningful for averaged cross spectra (m > 1)
- Single spectra have coherence ≈ 1.0 due to no ensemble averaging
- Low coherence may indicate non-linear coupling or noise dominance
"""
function coherence(cs::CrossSpectrum{T}) where {T}

    if is_single(cs)
        # For single cross spectra, coherence is theoretically 1.0 everywhere
        # but we can compute it anyway for diagnostic purposes
        coherence_values = Vector{T}(undef, length(cs.freq))

        for i in eachindex(cs.freq)
            cross_power_mag_sq = abs2(cs.power[i])
            coherence_values[i] = cross_power_mag_sq / (cs.ps1[i] * cs.ps2[i])
            # Clamp to [0, 1] to handle numerical issues
            coherence_values[i] = min(max(coherence_values[i], 0.0), 1.0)
        end

        return coherence_values

    else
        # For averaged cross spectra, coherence is meaningful
        coherence_values = Vector{T}(undef, length(cs.freq))

        for i in eachindex(cs.freq)
            cross_power_mag_sq = abs2(cs.power[i])
            coherence_values[i] = cross_power_mag_sq / (cs.ps1[i] * cs.ps2[i])
        end

        # For averaged spectra, apply bias correction
        # The expected coherence for pure noise is approximately 1/m
        # where m is the number of segments averaged
        bias_correction = 1.0 / cs.m

        # Subtract bias and ensure values stay in [0, 1]
        for i in eachindex(coherence_values)
            coherence_values[i] = max(coherence_values[i] - bias_correction, 0.0)
            coherence_values[i] = min(coherence_values[i], 1.0)
        end

        return coherence_values
    end
end

"""
    phase_lag(cs::CrossSpectrum) -> Vector

Calculate phase lag from cross spectrum.

Computes the phase difference between the two signals as a function of frequency.
Positive phase indicates the second signal lags the first.

# Arguments
- `cs::CrossSpectrum`: The cross spectrum to analyze

# Returns
- `Vector`: Phase lag in radians for each frequency bin

# Mathematical Details
- Phase lag = arg(P_xy(f)) = arctan(Im(P_xy)/Re(P_xy))
- Range: [-π, π] radians
- Positive values: signal 2 lags signal 1
- Negative values: signal 2 leads signal 1

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
phase = phase_lag(cs)
plot(cs.freq, phase, xlabel="Frequency (Hz)", ylabel="Phase Lag (rad)")
```
"""
function phase_lag(cs::CrossSpectrum)
    return angle.(cs.power)
end

"""
    time_lag(cs::CrossSpectrum) -> Vector

Calculate time lag from cross spectrum phase information.

Converts phase lag to time lag by dividing by 2πf. Positive values indicate
that the second signal lags behind the first signal.

# Arguments
- `cs::CrossSpectrum`: The cross spectrum to analyze

# Returns
- `Vector`: Time lag in seconds for each frequency bin

# Mathematical Details
- Time lag = phase_lag(f) / (2π × f)
- Units: seconds
- Positive values: signal 2 lags signal 1
- Negative values: signal 2 leads signal 1

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
tlag = time_lag(cs)
plot(cs.freq, tlag, xlabel="Frequency (Hz)", ylabel="Time Lag (s)")
```

# Notes
- Time lag becomes very large (and less meaningful) at low frequencies
- Consider filtering or masking very low frequency bins
- Assumes linear relationship between signals
"""
function time_lag(cs::CrossSpectrum)
    return phase_lag(cs) ./ (2π .* cs.freq)
end

"""
    noise_properties(cs::CrossSpectrum{T}) where T -> Dict{String, Any}

Comprehensive noise analysis for cross spectrum.

This function analyzes various noise properties of the cross spectrum including
theoretical and empirical noise levels, signal-to-noise statistics, and detection rates.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `Dict{String, Any}`: Dictionary containing noise analysis results with keys:
  - `"theoretical_noise"`: Theoretical Poisson noise level
  - `"noise_level_1"`, `"noise_level_2"`: Individual signal noise levels
  - `"mean_power"`, `"std_power"`: Power spectrum statistics
  - `"mean_snr"`: Average signal-to-noise ratio
  - `"significant_detections"`: Number of bins with S/N > 3
  - `"total_frequencies"`: Total number of frequency bins
  - `"high_freq_power"`: Mean power in high frequency region
  - `"segments_averaged"`: Number of segments averaged (if applicable)

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
props = noise_properties(cs)
println("Mean S/N: ", props["mean_snr"])
println("Significant detections: ", props["significant_detections"])
```
"""
function noise_properties(cs::CrossSpectrum{T}) where {T}

    # Individual noise levels
    if is_averaged(cs)
        noise1 = poisson_level(cs.norm; meanrate = cs.mean_rate1, n_ph = cs.nphots1 / cs.m)
        noise2 = poisson_level(cs.norm; meanrate = cs.mean_rate2, n_ph = cs.nphots2 / cs.m)
    else
        noise1 = poisson_level(cs.norm; meanrate = cs.mean_rate1, n_ph = cs.nphots1)
        noise2 = poisson_level(cs.norm; meanrate = cs.mean_rate2, n_ph = cs.nphots2)
    end

    cross_noise = theoretical_noise_level(cs)

    # Statistics
    mean_power = mean(abs.(cs.power))
    std_power = std(abs.(cs.power))

    # S/N analysis
    snr = signal_to_noise_ratio(cs)
    mean_snr = mean(snr)
    significant_bins = sum(snr .> 3.0)

    # High frequency noise estimate
    high_freq_idx = cs.freq .> (maximum(cs.freq) * 0.8)
    high_freq_power = length(high_freq_idx) > 0 ? mean(abs.(cs.power[high_freq_idx])) : NaN

    # Build properties dictionary
    properties = Dict(
        "theoretical_noise" => cross_noise,
        "noise_level_1" => noise1,
        "noise_level_2" => noise2,
        "mean_power" => mean_power,
        "std_power" => std_power,
        "mean_snr" => mean_snr,
        "significant_detections" => significant_bins,
        "total_frequencies" => length(cs.freq),
        "high_freq_power" => high_freq_power,
        "noise_to_signal_ratio" => cross_noise / mean_power,
        "mean_rate_1" => cs.mean_rate1,
        "mean_rate_2" => cs.mean_rate2,
        "is_averaged" => is_averaged(cs),
    )
    # Add averaging-specific properties
    if is_averaged(cs)
        properties["segments_averaged"] = cs.m
        properties["averaging_improvement"] = sqrt(cs.m)
    else
        properties["segments_averaged"] = 1
        properties["averaging_improvement"] = 1.0
    end

    return properties
end

"""
    significant_frequencies(cs::CrossSpectrum{T}, threshold::Real=3.0) where T -> Vector{T}

Return frequencies where signal detection is significant.

Identifies frequency bins where the signal-to-noise ratio exceeds a specified
threshold, indicating statistically significant signal detection.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze
- `threshold::Real=3.0`: S/N ratio threshold for significance

# Returns
- `Vector{T}`: Frequencies with significant signal detection

# Mathematical Details
- Uses `signal_to_noise_ratio()` for S/N calculation
- Default threshold of 3.0 corresponds to ~99.7% confidence
- Higher thresholds provide more conservative detection

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
sig_freqs = significant_frequencies(cs, 5.0)  # Very conservative
println("Significant frequencies: ", length(sig_freqs), " out of ", length(cs.freq))
```

# Applications
- Identifying QPO frequencies
- Filtering noise-dominated frequency bins
- Statistical significance testing
"""
function significant_frequencies(cs::CrossSpectrum{T}, threshold::Real = 3.0) where {T}
    snr = signal_to_noise_ratio(cs)
    significant_mask = snr .> threshold
    return cs.freq[significant_mask]
end

"""
    get_noise_level(cs::CrossSpectrum{T}, method::Symbol=:theoretical) where T -> T

Get noise level using different estimation methods.

Provides flexible access to different noise level estimation approaches,
either theoretical (based on Poisson statistics) or empirical (from high frequencies).

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze
- `method::Symbol=:theoretical`: Estimation method (`:theoretical` or `:empirical`)

# Returns
- `T`: Estimated noise level

# Methods
- `:theoretical`: Uses Poisson statistics and known count rates
- `:empirical`: Estimates from high frequency power levels

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
theoretical = get_noise_level(cs, :theoretical)
empirical = get_noise_level(cs, :empirical)
println("Theoretical: ", theoretical, ", Empirical: ", empirical)
```

# Throws
- `ArgumentError`: If unknown method specified
"""
function get_noise_level(cs::CrossSpectrum{T}, method::Symbol = :theoretical) where {T}
    if method == :theoretical
        return theoretical_noise_level(cs)
    elseif method == :empirical
        return white_noise_level(cs)
    else
        throw(ArgumentError("Unknown method: $method. Use :theoretical or :empirical"))
    end
end

"""
    quality_metrics(cs::CrossSpectrum{T}) where T -> Dict{String, Any}

Calculate quality metrics for cross spectrum analysis.

Computes various metrics to assess the quality and reliability of the cross
spectrum, including signal-to-noise statistics and dynamic range measures.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to analyze

# Returns
- `Dict{String, Any}`: Dictionary containing quality metrics:
  - `"mean_snr"`: Average signal-to-noise ratio
  - `"max_snr"`: Maximum signal-to-noise ratio
  - `"significant_fraction"`: Fraction of bins with S/N > 3
  - `"noise_level"`: Theoretical noise level
  - `"dynamic_range"`: Ratio of max to min power
  - `"is_averaged"`: Whether spectrum is averaged
  - `"segments_averaged"`: Number of averaged segments (if applicable)
  - `"expected_improvement"`: Expected S/N improvement from averaging

# Examples
```julia
cs = AveragedCrossSpectrum(lc1, lc2, 100.0)
metrics = quality_metrics(cs)
```

# Quality Indicators
- High `mean_snr`: Good overall signal quality
- High `significant_fraction`: Many reliable frequency bins
- High `dynamic_range`: Wide range of signal strengths
- `expected_improvement` ≈ √m for properly averaged spectra
"""
function quality_metrics(cs::CrossSpectrum{T}) where {T}
    snr = signal_to_noise_ratio(cs)

    metrics = Dict(
        "mean_snr" => mean(snr),
        "max_snr" => maximum(snr),
        "significant_fraction" => sum(snr .> 3.0) / length(snr),
        "noise_level" => theoretical_noise_level(cs),
        "dynamic_range" => maximum(abs.(cs.power)) / minimum(abs.(cs.power)),
        "is_averaged" => is_averaged(cs),
    )

    if is_averaged(cs)
        metrics["segments_averaged"] = cs.m
        metrics["expected_improvement"] = sqrt(cs.m)
    end

    return metrics
end

"""
    rebin(cs::CrossSpectrum{T}, df_new::Real) where T -> CrossSpectrum{T}

Rebin cross spectrum to a new frequency resolution.

Combines adjacent frequency bins to achieve lower frequency resolution and
improved signal-to-noise ratio through increased effective averaging.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to rebin
- `df_new::Real`: New frequency resolution (must be ≥ current resolution)

# Returns
- `CrossSpectrum{T}`: Rebinned cross spectrum with new frequency resolution

# Mathematical Details
- Rebin factor = round(df_new / current_df)
- New values are averages of original bins within each new bin
- Error propagation: σ_new = √(Σσ²) / n_bins for each new bin
- Updates k values to track effective number of averaged bins

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)  # df = 0.1 Hz
cs_rebinned = rebin(cs, 1.0)  # Rebin to 1.0 Hz resolution
println("Original bins: ", length(cs.freq))
println("Rebinned bins: ", length(cs_rebinned.freq))
```

# Notes
- Cannot increase frequency resolution (df_new < current_df)
- Some bins at the end may be dropped if not evenly divisible
- Preserves all spectrum metadata and normalization
"""
function rebin(cs::CrossSpectrum{T}, df_new::Real) where {T}
    if df_new < cs.df
        error("New frequency resolution must be >= current resolution")
    end

    rebin_factor = round(Int, df_new / cs.df)

    if rebin_factor == 1
        return cs
    end

    n_new = div(length(cs.freq), rebin_factor)
    if n_new == 0
        error("Rebinning factor too large for available frequencies")
    end

    freq_new = Vector{T}(undef, n_new)
    power_new = Vector{Complex{T}}(undef, n_new)
    ps1_new = Vector{T}(undef, n_new)
    ps2_new = Vector{T}(undef, n_new)
    power_err_new = isnothing(cs.power_err) ? nothing : Vector{T}(undef, n_new)

    for i = 1:n_new
        start_idx = (i - 1) * rebin_factor + 1
        end_idx = min(i * rebin_factor, length(cs.freq))

        freq_new[i] = mean(cs.freq[start_idx:end_idx])
        power_new[i] = mean(cs.power[start_idx:end_idx])
        ps1_new[i] = mean(cs.ps1[start_idx:end_idx])
        ps2_new[i] = mean(cs.ps2[start_idx:end_idx])

        if !isnothing(cs.power_err)
            power_err_new[i] =
                sqrt(sum(cs.power_err[start_idx:end_idx] .^ 2)) / (end_idx - start_idx + 1)
        end
    end

    k_new = if is_averaged(cs)
        if isa(cs.k, Int)
            cs.k * rebin_factor
        else
            rebin_factor
        end
    else
        rebin_factor
    end

    return CrossSpectrum{T}(
        freq_new,
        power_new,
        power_err_new,
        ps1_new,
        ps2_new,
        cs.norm,
        cs.power_type,
        df_new,
        cs.nphots1,
        cs.nphots2,
        cs.m,
        length(freq_new),
        k_new,
        cs.metadata1,
        cs.metadata2,
        cs.fullspec,
        cs.channels_overlap,
        cs.segment_size,
        cs.mean_rate1,
        cs.mean_rate2,
    )
end

"""
    rebin_log(cs::CrossSpectrum{T}; f::Real=0.01) where T -> CrossSpectrum{T}

Rebin cross spectrum using logarithmic frequency spacing.

Creates bins with logarithmically increasing widths, providing higher resolution
at low frequencies and lower resolution at high frequencies. This is particularly
useful for analyzing power-law spectra where features span multiple decades.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to rebin
- `f::Real=0.01`: Fractional frequency resolution (Δf/f), must be between 0 and 1

# Returns
- `CrossSpectrum{T}`: Logarithmically rebinned cross spectrum

# Mathematical Details
- Each bin spans frequency range [f_i, f_i(1+f)]
- Bin centers are geometric means: f_center = √(f_low × f_high)
- Number of bins ≈ log₁₀(f_max/f_min) / log₁₀(1+f)
- Preserves frequency-integrated power

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
cs_log = rebin_log(cs, f=0.1)  # 10% fractional resolution
plot(cs_log.freq, abs.(cs_log.power), xscale=:log10)
```

# Notes
- Automatically skips zero frequency if present
- Smaller f values give finer resolution but more bins
- Ideal for power-law noise analysis and broad-band features
"""
function rebin_log(cs::CrossSpectrum{T}; f::Real = 0.01) where {T}
    if f <= 0 || f >= 1
        error("Fractional frequency resolution f must be between 0 and 1")
    end

    start_idx = cs.freq[1] == 0 ? 2 : 1

    if start_idx >= length(cs.freq)
        error("Not enough frequency points for logarithmic rebinning")
    end

    freq_min = cs.freq[start_idx]
    freq_max = cs.freq[end]

    log_freq_min = log10(freq_min)
    log_freq_max = log10(freq_max)

    n_bins = floor(Int, (log_freq_max - log_freq_min) / log10(1 + f)) + 1

    if n_bins <= 1
        error("Not enough frequency range for logarithmic rebinning with f=$f")
    end

    freq_new = Vector{T}(undef, n_bins)
    power_new = Vector{Complex{T}}(undef, n_bins)
    ps1_new = Vector{T}(undef, n_bins)
    ps2_new = Vector{T}(undef, n_bins)
    power_err_new = isnothing(cs.power_err) ? nothing : Vector{T}(undef, n_bins)
    k_new = Vector{Int}(undef, n_bins)

    current_freq = freq_min

    for i = 1:n_bins
        freq_low = current_freq
        freq_high = current_freq * (1 + f)

        mask = (cs.freq .>= freq_low) .& (cs.freq .< freq_high)

        if i == n_bins
            mask = (cs.freq .>= freq_low) .& (cs.freq .<= freq_max)
        end

        indices = findall(mask)

        if isempty(indices)
            closest_idx = argmin(abs.(cs.freq .- (freq_low + freq_high) / 2))
            indices = [closest_idx]
        end

        freq_new[i] = sqrt(freq_low * freq_high)
        power_new[i] = mean(cs.power[indices])
        ps1_new[i] = mean(cs.ps1[indices])
        ps2_new[i] = mean(cs.ps2[indices])

        if is_averaged(cs)
            if isa(cs.k, Int)
                k_new[i] = cs.k * length(indices)
            else
                k_new[i] = sum(cs.k[indices])
            end
        else
            k_new[i] = length(indices)
        end

        if !isnothing(cs.power_err)
            power_err_new[i] = sqrt(sum(cs.power_err[indices] .^ 2)) / length(indices)
        end

        current_freq = freq_high
    end

    df_new = freq_new[2] - freq_new[1]

    return CrossSpectrum{T}(
        freq_new,
        power_new,
        power_err_new,
        ps1_new,
        ps2_new,
        cs.norm,
        cs.power_type,
        df_new,
        cs.nphots1,
        cs.nphots2,
        cs.m,
        length(freq_new),
        k_new,
        cs.metadata1,
        cs.metadata2,
        cs.fullspec,
        cs.channels_overlap,
        cs.segment_size,
        cs.mean_rate1,
        cs.mean_rate2,
    )
end

"""
    rebin(cs::CrossSpectrum{T}, rebin_factor::Int) where T -> CrossSpectrum{T}

Rebin cross spectrum by an integer factor.

Convenience method that combines adjacent frequency bins by the specified integer
factor, equivalent to calling `rebin(cs, cs.df * rebin_factor)`.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to rebin
- `rebin_factor::Int`: Integer rebinning factor (≥ 1)

# Returns
- `CrossSpectrum{T}`: Rebinned cross spectrum

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
cs_2x = rebin(cs, 2)  # Combine every 2 bins
cs_10x = rebin(cs, 10)  # Combine every 10 bins
```

# Notes
- Returns original spectrum unchanged if `rebin_factor ≤ 1`
- More intuitive than specifying exact frequency resolution
"""
function rebin(cs::CrossSpectrum{T}, rebin_factor::Int) where {T}
    if rebin_factor <= 1
        return cs
    end

    df_new = cs.df * rebin_factor
    return rebin(cs, df_new)
end

"""
    geometric_rebin(cs::CrossSpectrum{T}, factor::Real) where T -> CrossSpectrum{T}

Rebin cross spectrum with geometrically increasing bin widths.

Creates frequency bins that grow in width by a constant multiplicative factor,
providing a compromise between linear and logarithmic rebinning. Each successive
bin is `factor` times wider than the previous one.

# Arguments
- `cs::CrossSpectrum{T}`: The cross spectrum to rebin
- `factor::Real`: Geometric growth factor for bin widths (must be > 1.0)

# Returns
- `CrossSpectrum{T}`: Geometrically rebinned cross spectrum

# Mathematical Details
- First bin width = original df
- Subsequent bin widths: Δf_n = Δf_1 × factor^(n-1)
- Bin centers are arithmetic means of bin edges
- Preserves total integrated power

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
cs_geom = geometric_rebin(cs, 1.5)  # Each bin 50% wider than previous
cs_geom2 = geometric_rebin(cs, 2.0)  # Each bin twice as wide
```

# Applications
- Intermediate between linear and logarithmic spacing
- Useful for spectra with features at multiple scales
- Good for quasi-periodic oscillations (QPO) analysis

# Notes
- Factor must be > 1.0
- Automatically skips zero frequency if present
- Final bin may extend to maximum frequency
"""
function geometric_rebin(cs::CrossSpectrum{T}, factor::Real) where {T}
    if factor <= 1.0
        error("Geometric factor must be > 1.0")
    end

    start_idx = cs.freq[1] == 0 ? 2 : 1

    if start_idx >= length(cs.freq)
        error("Not enough frequency points for geometric rebinning")
    end

    freq_edges = [cs.freq[start_idx]]
    current_width = cs.df

    while freq_edges[end] < cs.freq[end]
        next_edge = freq_edges[end] + current_width
        if next_edge > cs.freq[end]
            freq_edges = push!(freq_edges, cs.freq[end])
            break
        end
        freq_edges = push!(freq_edges, next_edge)
        current_width *= factor
    end

    n_bins = length(freq_edges) - 1

    if n_bins <= 1
        error("Not enough frequency range for geometric rebinning")
    end

    freq_new = Vector{T}(undef, n_bins)
    power_new = Vector{Complex{T}}(undef, n_bins)
    ps1_new = Vector{T}(undef, n_bins)
    ps2_new = Vector{T}(undef, n_bins)
    power_err_new = isnothing(cs.power_err) ? nothing : Vector{T}(undef, n_bins)
    k_new = Vector{Int}(undef, n_bins)

    for i = 1:n_bins
        mask = (cs.freq .>= freq_edges[i]) .& (cs.freq .< freq_edges[i+1])

        if i == n_bins
            mask = (cs.freq .>= freq_edges[i]) .& (cs.freq .<= freq_edges[i+1])
        end

        indices = findall(mask)

        if isempty(indices)
            closest_idx = argmin(abs.(cs.freq .- (freq_edges[i] + freq_edges[i+1]) / 2))
            indices = [closest_idx]
        end

        freq_new[i] = (freq_edges[i] + freq_edges[i+1]) / 2
        power_new[i] = mean(cs.power[indices])
        ps1_new[i] = mean(cs.ps1[indices])
        ps2_new[i] = mean(cs.ps2[indices])

        if is_averaged(cs)
            if isa(cs.k, Int)
                k_new[i] = cs.k * length(indices)
            else
                k_new[i] = sum(cs.k[indices])
            end
        else
            k_new[i] = length(indices)
        end

        if !isnothing(cs.power_err)
            power_err_new[i] = sqrt(sum(cs.power_err[indices] .^ 2)) / length(indices)
        end
    end

    df_new = freq_new[2] - freq_new[1]

    return CrossSpectrum{T}(
        freq_new,
        power_new,
        power_err_new,
        ps1_new,
        ps2_new,
        cs.norm,
        cs.power_type,
        df_new,
        cs.nphots1,
        cs.nphots2,
        cs.m,
        length(freq_new),
        k_new,
        cs.metadata1,
        cs.metadata2,
        cs.fullspec,
        cs.channels_overlap,
        cs.segment_size,
        cs.mean_rate1,
        cs.mean_rate2,
    )
end

"""
    is_rebinned(cs::CrossSpectrum) -> Bool

Check whether a cross spectrum has been rebinned.

Determines if the cross spectrum has undergone any rebinning operation by
examining the `k` parameter, which tracks the number of original frequency
bins combined into each current bin.

# Arguments
- `cs::CrossSpectrum`: Cross spectrum to check

# Returns
- `Bool`: `true` if rebinned, `false` if original resolution

# Mathematical Details
- Original spectrum: k = 1 (scalar)
- Uniformly rebinned: k > 1 (scalar)  
- Non-uniformly rebinned: k is Vector with varying values

# Examples
```julia
cs_orig = CrossSpectrum(lc1, lc2)
println(is_rebinned(cs_orig))  # false

cs_rebin = rebin(cs_orig, 5)
println(is_rebinned(cs_rebin))  # true

cs_log = rebin_log(cs_orig)
println(is_rebinned(cs_log))  # true
```

# Applications
- Quality control in analysis pipelines
- Determining appropriate statistical methods
- Choosing correct error propagation formulas
"""
function is_rebinned(cs::CrossSpectrum)
    return isa(cs.k, Vector) || (isa(cs.k, Int) && cs.k > 1)
end

"""
    effective_samples_per_bin(cs::CrossSpectrum) -> Union{Int, Vector{Int}}

Calculate effective number of independent samples per frequency bin.

Computes the total number of independent measurements contributing to each
frequency bin, accounting for both segment averaging and frequency rebinning.
This is crucial for proper statistical analysis and error estimation.

# Arguments
- `cs::CrossSpectrum`: Cross spectrum to analyze

# Returns
- `Int`: Effective samples per bin (if uniform)
- `Vector{Int}`: Effective samples for each bin (if non-uniform rebinning)

# Mathematical Details
- For averaged spectrum: N_eff = m × k
- For single spectrum: N_eff = k
- Where m = number of segments, k = rebinning factor per bin

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)  # Single segment
n_eff = effective_samples_per_bin(cs)  # Returns 1

cs_avg = average([cs1, cs2, cs3])  # 3 segments
n_eff_avg = effective_samples_per_bin(cs_avg)  # Returns 3

cs_rebin = rebin(cs_avg, 5)  # Rebin by factor 5
n_eff_rebin = effective_samples_per_bin(cs_rebin)  # Returns 15
```

# Applications
- Chi-squared goodness-of-fit testing
- Confidence interval calculation
- Model parameter uncertainty estimation
- Detection significance assessment

# Notes
- Higher values indicate better statistical precision
- Vector return indicates non-uniform rebinning (e.g., logarithmic)
- Essential for proper statistical interpretation
"""
function effective_samples_per_bin(cs::CrossSpectrum)
    if is_averaged(cs)
        if isa(cs.k, Int)
            return cs.m * cs.k
        else
            return cs.m .* cs.k
        end
    else
        if isa(cs.k, Int)
            return cs.k
        else
            return cs.k
        end
    end
end

"""
    adaptive_rebin(cs::CrossSpectrum{T}, target_snr::Real=3.0, max_rebin_factor::Int=10) where T -> CrossSpectrum{T}

Automatically rebin cross spectrum to achieve target signal-to-noise ratio.

Iteratively increases the rebinning factor until the desired SNR is achieved
or the maximum rebinning factor is reached. This is useful for optimizing
the trade-off between frequency resolution and statistical significance.

# Arguments
- `cs::CrossSpectrum{T}`: Cross spectrum to adaptively rebin
- `target_snr::Real=3.0`: Target signal-to-noise ratio
- `max_rebin_factor::Int=10`: Maximum allowed rebinning factor

# Returns
- `CrossSpectrum{T}`: Rebinned cross spectrum meeting SNR requirement

# Algorithm
1. Calculate current SNR for each frequency bin
2. If mean SNR ≥ target_snr, return original spectrum
3. Otherwise, try rebinning factors 2, 3, ..., max_rebin_factor
4. Return first rebinning that achieves target SNR
5. If none succeed, return spectrum with maximum rebinning

# Examples
```julia
cs_noisy = CrossSpectrum(lc1, lc2)  # Low SNR spectrum
cs_clean = adaptive_rebin(cs_noisy, target_snr=5.0)

# For detection work
cs_detect = adaptive_rebin(cs, target_snr=3.0)  # 3σ detection
```

# Applications
- Automatic optimization of frequency resolution vs. sensitivity
- Standardizing analysis pipelines for different data qualities
- Preparing spectra for feature detection algorithms
- Quality-driven data reduction

# Notes
- Uses linear rebinning (uniform frequency spacing)
- SNR calculation depends on available error estimates
- May significantly reduce frequency resolution
- Consider scientific requirements before applying
"""
function adaptive_rebin(
    cs::CrossSpectrum{T},
    target_snr::Real = 3.0,
    max_rebin_factor::Int = 10,
) where {T}
    current_snr = signal_to_noise_ratio(cs)

    needs_rebinning = current_snr .< target_snr

    if !any(needs_rebinning)
        return cs
    end

    for factor = 2:max_rebin_factor
        rebinned_cs = rebin(cs, factor)
        new_snr = signal_to_noise_ratio(rebinned_cs)

        if mean(new_snr) >= target_snr
            return rebinned_cs
        end
    end

    return rebin(cs, max_rebin_factor)
end

"""
    Base.show(io::IO, ::MIME"text/plain", cs::CrossSpectrum)

Display cross spectrum information in a readable format.

Provides a comprehensive summary of the cross spectrum properties,
including frequency information, averaging details, and analysis parameters.
The display format differs based on whether the spectrum is averaged or single.

# Arguments
- `io::IO`: Output stream
- `::MIME"text/plain"`: MIME type for plain text display
- `cs::CrossSpectrum`: Cross spectrum to display

# Display Information

## For Averaged Spectra
- Type and element type
- Number of frequency bins
- Number of segments averaged
- Segment size and mean count rates
- Normalization method

## For Single Spectra  
- Type and element type
- Number of frequency bins
- Normalization method

# Examples
```julia
cs = CrossSpectrum(lc1, lc2)
println(cs)
# Output:
# CrossSpectrum{Float64} (Single)
#   Frequencies: 1024
#   Normalization: leahy

cs_avg = average([cs1, cs2, cs3])
println(cs_avg)
# Output:
# CrossSpectrum{Float64} (Averaged)
#   Frequencies: 1024
#   Segments averaged: 3
#   Segment size: 2048
#   Mean rates: 150.0, 200.0
#   Normalization: leahy
```

# Notes
- Automatically called by `println()` and REPL display
- Helps distinguish between single and averaged spectra
- Provides key information for analysis interpretation
"""
function Base.show(io::IO, ::MIME"text/plain", cs::CrossSpectrum)
    if is_averaged(cs)
        println(io, "CrossSpectrum{$(eltype(cs.freq))} (Averaged)")
        println(io, "  Frequencies: $(length(cs.freq))")
        println(io, "  Segments averaged: $(cs.m)")
        println(io, "  Segment size: $(cs.segment_size)")
        println(io, "  Mean rates: $(cs.mean_rate1), $(cs.mean_rate2)")
        println(io, "  Normalization: $(cs.norm)")
    else
        println(io, "CrossSpectrum{$(eltype(cs.freq))} (Single)")
        println(io, "  Frequencies: $(length(cs.freq))")
        println(io, "  Normalization: $(cs.norm)")
    end
end