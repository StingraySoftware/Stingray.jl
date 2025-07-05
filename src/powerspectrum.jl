"""
Abstract type representing a power spectrum, which characterizes the distribution 
of power across different frequencies in a signal.

Subtypes include:
- PowerSpectrum{T}: Represents a power spectrum for a single signal segment
- AveragedPowerspectrum{T}: Represents a power spectrum averaged over multiple segments

# Type Parameters
- `T`: The numeric type for frequency and power values (typically Float64)
"""
abstract type AbstractPowerSpectrum{T} end

"""
    PowerSpectrum{T} <: AbstractPowerSpectrum{T}

Represents a power spectrum computed from a single signal segment.

# Fields
- `freqs::Vector{T}`: Vector of frequency values 
- `power::Vector{T}`: Corresponding power values for each frequency
- `power_errors::Vector{T}`: Uncertainty estimates for power at each frequency
- `dt::T`: Time sampling interval
- `n::Int`: Number of data points in the original segment

# Examples
```julia
# Create a power spectrum from raw data
freqs = [0.1, 0.2, 0.3, 0.4, 0.5]
power = [1.2, 2.3, 1.8, 0.9, 1.1]
errors = [0.1, 0.2, 0.15, 0.08, 0.09]
ps = PowerSpectrum(freqs, power, errors, 0.001, 1024)
```
"""
struct PowerSpectrum{T} <: AbstractPowerSpectrum{T}
    freqs::Vector{T}
    power::Vector{T}
    power_errors::Vector{T}
    dt::T
    n::Int
end

"""
    AveragedPowerspectrum{T} <: AbstractPowerSpectrum{T}

Represents a power spectrum computed by averaging multiple signal segments using the 
Bartlett method for reduced variance.

# Fields
- `freqs::Vector{T}`: Vector of frequency values
- `power::Vector{T}`: Average power values for each frequency
- `power_errors::Vector{T}`: Standard error of power estimates
- `m::Int`: Number of segments used in averaging
- `n::Int`: Number of data points in each segment
- `dt::T`: Time sampling interval

# Examples
```julia
# Create an averaged power spectrum
freqs = collect(0.1:0.1:5.0)
power = rand(length(freqs))
errors = rand(length(freqs)) * 0.1
aps = AveragedPowerspectrum(freqs, power, errors, 32, 1024, 0.001)
```
"""
struct AveragedPowerspectrum{T} <: AbstractPowerSpectrum{T}
    freqs::Vector{T}
    power::Vector{T} 
    power_errors::Vector{T}
    m::Int
    n::Int
    dt::T
end

"""
    powerspectrum(events::EventList; dt::Real=0.001, segment_size::Int=length(events.times) ÷ 32, norm::Symbol=:frac) -> AveragedPowerspectrum

Compute the power spectrum for an EventList using segmented averaging (Bartlett method).

The function converts the event times into a binned light curve and then computes the 
power spectrum using FFT with Hanning windowing to reduce spectral leakage.

# Arguments
- `events::EventList`: Input event list containing times and optional energies
- `dt::Real=0.001`: Time sampling interval for binning (default: 0.001)
- `segment_size::Int=length(events.times) ÷ 32`: Size of each segment for FFT computation
- `norm::Symbol=:frac`: Normalization method for power spectrum
    - `:leahy`: Leahy normalization (multiplies power by 2)
    - `:frac`: Fractional RMS normalization (normalizes by mean power)
    - `:abs`: Absolute normalization (no normalization applied)

# Returns
- `AveragedPowerspectrum{T}`: Power spectrum computed from segmented data where T matches the type of events.times

# Throws
- `ArgumentError`: If segment_size ≤ 0 or if there are insufficient data points

# Examples
```julia
# Basic usage with default parameters
events = EventList(times, energies)
ps = powerspectrum(events)

# Custom parameters
ps = powerspectrum(events, dt=0.01, segment_size=1024, norm=:leahy)

# High-resolution analysis
ps = powerspectrum(events, dt=0.0001, segment_size=4096, norm=:frac)
```

# Notes
- The function automatically removes the DC component (frequency = 0) from the output
- Hanning windowing is applied to reduce spectral leakage
- Power errors are computed as standard error of the mean across segments
- Frequency resolution is determined by Δf = 1/(segment_size × dt)
"""
function powerspectrum(events::EventList; dt::Real = 0.001,
                      segment_size::Int = length(events.times) ÷ 32, 
                      norm::Symbol = :frac)
    T = eltype(events.times)
    
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive, got $segment_size"))
    end
    if length(events.times) < segment_size
        throw(ArgumentError("Not enough data points ($(length(events.times))) for the specified segment size ($segment_size)"))
    end
    
    return powerspectrum_averaged(events, convert(T, dt), segment_size, norm)
end

"""
    powerspectrum(lc::LightCurve; segment_size::Int=length(lc.time) ÷ 32, norm::Symbol=:frac) -> AveragedPowerspectrum

Compute the power spectrum for a LightCurve using segmented averaging (Bartlett method).

The function directly uses the binned light curve data and computes the power spectrum 
using FFT with Hanning windowing.

# Arguments
- `lc::LightCurve`: Input light curve with time and count data
- `segment_size::Int=length(lc.time) ÷ 32`: Size of each segment for FFT computation
- `norm::Symbol=:frac`: Normalization method (see powerspectrum(::EventList) for options)

# Returns
- `AveragedPowerspectrum{T}`: Power spectrum where T matches the type of lc.time

# Throws
- `ArgumentError`: If segment_size ≤ 0 or if there are insufficient data points

# Examples
```julia
# Basic usage
lc = LightCurve(times, counts, metadata)
ps = powerspectrum(lc)

# Custom segment size for better frequency resolution
ps = powerspectrum(lc, segment_size=2048, norm=:leahy)
```

# Notes
- Uses the bin_size from light curve metadata as the time sampling interval
- Automatically handles the conversion of count data to appropriate numeric type
"""
function powerspectrum(lc::LightCurve; segment_size::Int = length(lc.time) ÷ 32, 
                      norm::Symbol = :frac)
    T = eltype(lc.time)
    
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive, got $segment_size"))
    end
    if length(lc.time) < segment_size
        throw(ArgumentError("Not enough data points ($(length(lc.time))) for the specified segment size ($segment_size)"))
    end
    
    # Use the bin size from the light curve metadata as dt
    dt = T(lc.metadata.bin_size)
    
    return powerspectrum_averaged(lc, dt, segment_size, norm)
end

"""
    powerspectrum_averaged(events::EventList, dt::T, segment_size::Int, norm::Symbol) where T -> AveragedPowerspectrum{T}

Internal function to compute averaged power spectrum from EventList segments.

This function first converts the event list to a binned light curve, then computes
the power spectrum using the Bartlett method with Hanning windowing.

# Arguments
- `events::EventList`: Input event list
- `dt::T`: Time sampling interval
- `segment_size::Int`: Size of each segment for FFT
- `norm::Symbol`: Normalization method

# Returns
- `AveragedPowerspectrum{T}`: Computed power spectrum

# Algorithm for devlopers :)
1. Create time bins from minimum to maximum event time
2. Bin events into histogram counts
3. Apply segmented FFT analysis with Hanning windowing
4. Average power spectra across segments
5. Compute standard errors and apply normalization
"""
function powerspectrum_averaged(events::EventList, 
                               dt::T, 
                               segment_size::Int, 
                               norm::Symbol) where T
    
    if length(events.times) == 0
        throw(ArgumentError("EventList is empty"))
    end
    
    # Create time bins from minimum to maximum event time
    t_min = minimum(events.times)
    t_max = maximum(events.times)
    
    # Create bin edges
    n_bins = round(Int, (t_max - t_min) / dt) + 1
    bin_edges = range(t_min, step=dt, length=n_bins+1)
    
    # Bin the events to create a light curve
    counts = zeros(Int, n_bins)
    for t in events.times
        if t >= t_min && t <= t_max
            bin_idx = min(n_bins, max(1, round(Int, (t - t_min) / dt) + 1))
            counts[bin_idx] += 1
        end
    end
    
    # Now compute power spectrum from the binned data
    return compute_powerspectrum_from_counts(counts, dt, segment_size, norm)
end

"""
    powerspectrum_averaged(lc::LightCurve, dt::T, segment_size::Int, norm::Symbol) where T -> AveragedPowerspectrum{T}

Internal function to compute averaged power spectrum from LightCurve segments.

This function directly uses the light curve count data to compute the power spectrum
using the Bartlett method with Hanning windowing.

# Arguments
- `lc::LightCurve`: Input light curve
- `dt::T`: Time sampling interval
- `segment_size::Int`: Size of each segment for FFT
- `norm::Symbol`: Normalization method

# Returns
- `AveragedPowerspectrum{T}`: Computed power spectrum
"""
function powerspectrum_averaged(lc::LightCurve, 
                               dt::T, 
                               segment_size::Int, 
                               norm::Symbol) where T
    
    # Convert counts to float for FFT processing
    counts_float = convert(Vector{Float64}, lc.counts)
    
    return compute_powerspectrum_from_counts(counts_float, dt, segment_size, norm)
end

"""
    compute_powerspectrum_from_counts(counts::Vector{<:Real}, dt::T, segment_size::Int, norm::Symbol) where T -> AveragedPowerspectrum{T}

Core function to compute power spectrum from count data using the Bartlett method.

This is the main computational engine that implements the segmented FFT analysis
with Hanning windowing and statistical error estimation.

# Arguments
- `counts::Vector{<:Real}`: Time series count data
- `dt::T`: Time sampling interval
- `segment_size::Int`: Size of each segment for FFT
- `norm::Symbol`: Normalization method

# Returns
- `AveragedPowerspectrum{T}`: Computed power spectrum

# Algorithm for devlopers :)
1. Divide the time series into overlapping or non-overlapping segments
2. Apply Hanning window to each segment to reduce spectral leakage
3. Compute FFT of each windowed segment
4. Calculate power spectrum (|FFT|²) for each segment
5. Average power spectra across all segments
6. Compute standard errors as standard deviation / √(number of segments)
7. Apply specified normalization
8. Remove DC component (frequency = 0)

# Throws
- `ArgumentError`: If signal is too short for the specified segment size
"""
function compute_powerspectrum_from_counts(counts::Vector{<:Real}, dt::T, 
                                         segment_size::Int, norm::Symbol) where T
    
    n_segments = length(counts) ÷ segment_size
    if n_segments == 0
        throw(ArgumentError("Signal too short ($(length(counts)) points) for the specified segment size ($segment_size)"))
    end
    
    # Calculate mean count rate for proper fractional normalization
    mean_count_rate = mean(counts)
    
    # Create segments
    segments = Vector{Vector{Float64}}()
    for i in 1:n_segments
        start_idx = (i-1) * segment_size + 1
        end_idx = i * segment_size
        push!(segments, counts[start_idx:end_idx])
    end
    
    # Generate Hanning window
    window = hanning(segment_size)
    
    # Compute window normalization factor
    window_norm = sum(window.^2) / segment_size
    
    # Apply windowing and compute FFTs
    windowed_segments = Vector{Vector{ComplexF64}}()
    for segment in segments
        # Remove DC component before windowing
        segment_mean = mean(segment)
        detrended = segment .- segment_mean
        windowed = detrended .* window
        push!(windowed_segments, fft(windowed))
    end
    
    # Compute power spectra (only positive frequencies, excluding DC)
    nyquist_index = segment_size ÷ 2 + 1
    powers = Vector{Vector{Float64}}()
    
    for fft_result in windowed_segments
        # Skip DC component (index 1), take positive frequencies
        power_spectrum = abs2.(fft_result[2:nyquist_index])
        
        # Apply proper scaling for power spectrum
        # Factor of 2 accounts for negative frequencies (except DC and Nyquist)
        power_spectrum .*= 2.0 / (segment_size * window_norm)
        
        # Handle Nyquist frequency (if segment_size is even)
        if segment_size % 2 == 0
            power_spectrum[end] *= 0.5  # Nyquist frequency appears only once
        end
        
        push!(powers, power_spectrum)
    end
    
    # Average the power spectra
    n_freq_bins = length(powers[1])
    avg_power = zeros(Float64, n_freq_bins)
    for power_spectrum in powers
        avg_power .+= power_spectrum
    end
    avg_power ./= length(powers)
    
    # Compute frequency array (excluding DC component)
    freqs = collect(1:(n_freq_bins)) .* (1.0 / (segment_size * dt))
    
    # Calculate standard errors
    if length(powers) > 1
        power_matrix = hcat(powers...)
        power_stds = [std(power_matrix[i, :]) for i in 1:size(power_matrix, 1)]
        power_errors = power_stds ./ sqrt(length(segments))
    else
        # Single segment case - use sqrt(power) as approximation
        power_errors = sqrt.(avg_power) ./ sqrt(length(segments))
    end
    
    # Apply normalization with proper mean count rate
    normalized_power, normalized_errors = normalize_power_with_errors(avg_power, power_errors, norm, mean_count_rate)
    
    return AveragedPowerspectrum{T}(
        convert(Vector{T}, freqs),
        convert(Vector{T}, normalized_power),
        convert(Vector{T}, normalized_errors),
        length(segments),
        segment_size,
        dt
    )
end

"""
    hanning(N::Int) -> Vector{Float64}

Generate a Hanning (Hann) window function for signal tapering.

The Hanning window reduces spectral leakage by smoothly tapering the signal
to zero at the boundaries, which is crucial for accurate power spectral analysis.

# Arguments
- `N::Int`: Length of the window (number of points)

# Returns
- `Vector{Float64}`: Hanning window coefficients

# Formula
The Hanning window is defined as:
w(n) = 0.5 * (1 - cos(2π * n / (N-1))) for n = 0, 1, ..., N-1

# Examples
```julia
# Create a 1024-point Hanning window
window = hanning(1024)

# Apply to a signal
windowed_signal = signal .* window
```

# Notes
- Also known as the Hann window
- Provides good balance between main lobe width and side lobe suppression
- Commonly used in power spectral density estimation
"""
function hanning(N::Int)
    if N <= 0
        throw(ArgumentError("Window length must be positive, got $N"))
    end
    if N == 1
        return [1.0]
    end
    return [0.5 * (1 - cos(2π * n / (N-1))) for n in 0:N-1]
end

"""
    normalize_power(power::AbstractVector, norm::Symbol) -> Vector{Float64}

Normalize power spectrum based on specified method.

Different normalization methods are used depending on the analysis requirements
and the type of signals being studied.

# Arguments
- `power::AbstractVector`: Raw power spectrum values
- `norm::Symbol`: Normalization method

# Normalization Methods
- `:leahy`: Leahy normalization (power × 2) - preserves Poisson statistics
- `:frac`: Fractional RMS normalization (power / mean(power)) - emphasizes variability
- `:abs`: Absolute normalization (no change) - preserves original units

# Returns
- `Vector{Float64}`: Normalized power spectrum

# Examples
```julia
raw_power = [1.2, 2.3, 1.8, 0.9, 1.1]
leahy_power = normalize_power(raw_power, :leahy)
frac_power = normalize_power(raw_power, :frac)
abs_power = normalize_power(raw_power, :abs)
```

# Throws
- `ArgumentError`: If normalization method is not recognized
"""
function normalize_power_with_errors(power::AbstractVector, errors::AbstractVector, norm::Symbol, mean_count_rate::Real=1.0)
    if norm == :leahy
        # Leahy normalization: multiply by 2
        return power .* 2.0, errors .* 2.0
    elseif norm == :frac
        # Fractional RMS normalization: normalize by mean count rate
        # For white noise, this should give a mean power of 2/mean_count_rate
        # But we want the noise level to be ~1, so we need to scale appropriately
        
        # The correct fractional normalization is:
        # P_frac = P_raw * 2 / mean_count_rate
        # This gives noise level of 2/mean_count_rate for white noise
        # To get noise level of 1, we multiply by mean_count_rate/2
        
        # However, looking at the standard definition of fractional RMS:
        # We want the power to be normalized so that white noise has power ~1
        # This means: P_frac = P_raw * 2 / mean_count_rate
        
        if mean_count_rate > 0
            scale_factor = 2.0 / mean_count_rate
            normalized_power = power .* scale_factor
            normalized_errors = errors .* scale_factor
        else
            # Fallback if mean_count_rate is not provided or is zero
            normalized_power = power ./ 2.0
            normalized_errors = errors ./ 2.0
        end
        
        return normalized_power, normalized_errors
    elseif norm == :abs
        # Absolute normalization: no change
        return copy(power), copy(errors)
    else
        throw(ArgumentError("Unknown normalization method: $norm. Supported methods are :leahy, :frac, :abs"))
    end
end
function normalize_power(power::AbstractVector, norm::Symbol)
    normalized_power, _ = normalize_power_with_errors(power, zeros(length(power)), norm)
    return normalized_power
end
# Accessor functions
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
freqs(ps::AbstractPowerSpectrum) = ps.freqs

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
errors(ps::AbstractPowerSpectrum) = ps.power_errors

# Base methods
"""
    length(ps::AbstractPowerSpectrum) -> Int

Get the number of frequency bins in the power spectrum.
"""
Base.length(ps::AbstractPowerSpectrum) = length(ps.freqs)

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
Base.getindex(ps::AbstractPowerSpectrum, i) = (ps.freqs[i], ps.power[i], ps.power_errors[i])

# Show methods
"""
Display a compact representation of a PowerSpectrum.
"""
function Base.show(io::IO, ps::PowerSpectrum{T}) where T
    df = length(ps.freqs) > 1 ? ps.freqs[2] - ps.freqs[1] : 0.0
    print(io, "PowerSpectrum{$T}(n=$(ps.n), df=$(df))")
end

"""
Display a compact representation of an AveragedPowerspectrum.
"""
function Base.show(io::IO, ps::AveragedPowerspectrum{T}) where T
    df = length(ps.freqs) > 1 ? ps.freqs[2] - ps.freqs[1] : 0.0
    print(io, "AveragedPowerspectrum{$T}(n=$(ps.n), m=$(ps.m), df=$(df))")
end

"""
    validate(ps::AbstractPowerSpectrum) -> Bool

Validate the power spectrum structure for consistency and physical validity.

# Arguments
- `ps::AbstractPowerSpectrum`: Power spectrum to validate

# Returns
- `Bool`: true if valid (throws exception if invalid)

# Throws
- `ArgumentError`: If validation fails with descriptive error message

# Validation Checks
- Frequency, power, and error arrays have same length
- Power values are non-negative
- Error values are non-negative
- Frequency values are monotonically increasing

# Examples
```julia
ps = powerspectrum(events)
validate(ps)  # Returns true or throws error
```
"""
function validate(ps::AbstractPowerSpectrum)
    if length(ps.freqs) != length(ps.power)
        throw(ArgumentError("Frequency and power arrays must have the same length: $(length(ps.freqs)) vs $(length(ps.power))"))
    end
    if length(ps.power) != length(ps.power_errors)
        throw(ArgumentError("Power and error arrays must have the same length: $(length(ps.power)) vs $(length(ps.power_errors))"))
    end
    if any(x -> x < 0, ps.power)
        throw(ArgumentError("Power values must be non-negative"))
    end
    if any(x -> x < 0, ps.power_errors)
        throw(ArgumentError("Error values must be non-negative"))
    end
    if length(ps.freqs) > 1 && !issorted(ps.freqs)
        throw(ArgumentError("Frequency values must be monotonically increasing"))
    end
    return true
end
