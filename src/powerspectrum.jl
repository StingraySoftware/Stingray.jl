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
    bin_size = lc.metadata.bin_size
    n_bins_per_segment = round(Int, segment_size / bin_size)
    
    # Input validation
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
    
    @info "Processing light curve" total_bins=length(lc.counts) segment_size=segment_size bins_per_segment=n_bins_per_segment gti_shape=size(gtis)
    
    # Use the appropriate generator based on whether data is binned
    segment_generator = generate_indices_of_segment_boundaries_binned(
        lc.time, gtis, segment_size, dt=bin_size
    )
    
    # Initialize arrays
    freqs = fftfreq(n_bins_per_segment, 1/bin_size)
    pos_freq_idx = positive_fft_bins(n_bins_per_segment; include_zero=false)
    freqs = freqs[pos_freq_idx]
    df = freqs[2] - freqs[1]
    
    total_power = zeros(T, length(pos_freq_idx))
    total_counts = 0
    n_segments_used = 0
    
    # Process each segment using the generator
    for (start_time, stop_time, start_idx, stop_idx) in segment_generator
        # Verify segment length
        segment_length = stop_idx - start_idx
        if segment_length != n_bins_per_segment
            @debug "Skipping segment: wrong length" segment_length n_bins_per_segment
            continue
        end
        
        # Extract segment data
        segment_counts = @view lc.counts[start_idx+1:stop_idx]
        
        # Skip if zero counts
        segment_sum = sum(segment_counts)
        if segment_sum == 0
            @debug "Skipping segment: zero counts" start_time stop_time
            continue
        end
        
        # Calculate FFT and power
        ft = fft(segment_counts)
        unnorm_power = abs2.(ft[pos_freq_idx])
        
        # Normalize power using functions from Fourier.jl
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
    
    @info "Segment processing complete" segments_used=n_segments_used
    
    if n_segments_used == 0
        throw(ArgumentError("No valid segments found"))
    end
    
    # Calculate final results
    avg_power = total_power ./ n_segments_used
    mean_rate = total_counts / (n_segments_used * segment_size)
    
    # Calculate errors using functions from Fourier.jl
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