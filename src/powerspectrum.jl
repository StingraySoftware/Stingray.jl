"""
Abstract type representing a power spectrum, which characterizes the distribution 
of power across different frequencies in a signal.

Subtypes include:
- PowerSpectrum: Represents a power spectrum for a single signal segment
- AveragedPowerspectrum: Represents a power spectrum averaged over multiple segments
"""
abstract type AbstractPowerSpectrum end

"""
Represents a power spectrum computed from a single signal segment.

## Fields
- `freqs`: Vector of frequency values 
- `power`: Corresponding power values for each frequency
- `power_errors`: Uncertainty estimates for power at each frequency
- `dt`: Time sampling interval
- `n`: Number of data points in the original segment
"""
struct PowerSpectrum <: AbstractPowerSpectrum
    freqs::Vector{Float64}
    power::Vector{Float64}
    power_errors::Vector{Float64}
    dt::Float64
    n::Int
end

"""
Represents a power spectrum computed by averaging multiple signal segments.

## Fields
- `freqs`: Vector of frequency values
- `power`: Average power values for each frequency
- `power_errors`: Standard error of power estimates
- `m`: Number of segments used in averaging
- `n`: Number of data points in each segment
- `dt`: Time sampling interval
"""
struct AveragedPowerspectrum <: AbstractPowerSpectrum
    freqs::Vector{Float64}
    power::Vector{Float64} 
    power_errors::Vector{Float64}
    m::Int
    n::Int
    dt::Float64
end

"""
Compute the power spectrum for a list of events using segmented averaging.

## Arguments
- `events::EventList{T}`: Input event list containing times and optional energies
- `dt::Real=0.001`: Time sampling interval (default: 0.001)
- `segment_size::Int=length(events.times) ÷ 32`: Size of each segment for FFT (default: 1/32 of total data)
- `norm::String="frac"`: Normalization method for power spectrum
    - "leahy": Multiplies power by 2
    - "frac": Normalizes power by mean power
    - "abs": No normalization

## Returns
- `AveragedPowerspectrum`: Power spectrum computed from segmented data

## Throws
- `AssertionError` if segment size is invalid or insufficient data points

## Examples
```julia
events = EventList(times, energies)
ps = powerspectrum(events, dt=0.01, segment_size=1024, norm="leahy")
```
"""
function powerspectrum(events::EventList{T}; dt::Real = 0.001,
                      segment_size::Int = length(events.times) ÷ 32, 
                      norm::String = "frac") where T
    @assert segment_size > 0 "Segment size must be positive"
    @assert length(events.times) >= segment_size "Not enough data points"
    return _powerspectrum_averaged(events, float(dt), segment_size, norm)
end

"""
Compute averaged power spectrum from event list segments.

## Arguments
- `events::EventList`: Input event list
- `dt::Float64`: Time sampling interval
- `segment_size::Int`: Number of points in each segment
- `norm::String`: Normalization method for power spectrum

## Returns
- `AveragedPowerspectrum`: Power spectrum averaged over multiple segments

## Notes
- Prioritizes using event energies if available, otherwise uses event times
- Computes FFT for each segment and averages power
- Calculates frequency-dependent power with specified normalization
"""
function _powerspectrum_averaged(events::EventList, 
                               dt::Float64, 
                               segment_size::Int, 
                               norm::String)
    signal = events.times
    
    if !all(iszero, events.energies)
        signal = events.energies
    end
    
    n_segments = length(signal) ÷ segment_size
    segments = [signal[(i-1)*segment_size+1 : i*segment_size] for i in 1:n_segments]
    
    ffts = [fft(segment) for segment in segments]
    
    nyquist_index = segment_size ÷ 2 + 1
    powers = [abs2.(fft[1:nyquist_index]) for fft in ffts]
    
    avg_power = mean(powers)
    
    freqs = FFTW.rfftfreq(segment_size, 1/dt)
    @assert length(freqs) == length(avg_power) "Frequency and power arrays must match"
    
    normalized_power = _normalize_power(avg_power, norm)
    
    power_errors = std(powers) ./ sqrt(length(segments))
    
    positive_indices = 2:length(freqs)
    
    return AveragedPowerspectrum(
        freqs[positive_indices],
        normalized_power[positive_indices],
        power_errors[positive_indices],
        length(segments),
        segment_size,
        dt
    )
end

"""
Generate a Hanning window function for signal tapering.

## Arguments
- `N::Int`: Length of the window

## Returns
- `Vector{Float64}`: Hanning window weights

## Examples
```julia
window = hanning(1024)  # Creates a Hanning window of length 1024
```
"""
function hanning(N::Int)
    return [0.5 * (1 - cos(2π * n / (N-1))) for n in 0:N-1]
end

"""
Prepare a signal into overlapping or non-overlapping segments.

## Arguments
- `signal::Vector{Float64}`: Input signal
- `segment_size::Int`: Size of each segment

## Returns
- `Vector{Vector{Float64}}`: Signal divided into segments
"""
function _prepare_segments(signal::Vector{Float64}, segment_size::Int)
    n_segments = length(signal) ÷ segment_size
    segments = [signal[(i-1)*segment_size+1 : i*segment_size] for i in 1:n_segments]
    return segments
end

"""
Normalize power spectrum based on specified method.

## Arguments
- `power::Vector{Float64}`: Input power values
- `norm::String`: Normalization method
    - "leahy": Multiply by 2
    - "frac": Divide by mean power
    - "abs": No normalization

## Returns
- `Vector{Float64}`: Normalized power values

## Throws
- `ArgumentError` for unknown normalization method
"""
function _normalize_power(power::Vector{Float64}, norm::String)
    if norm == "leahy"
        return power .* 2
    elseif norm == "frac"
        return power ./ mean(power)
    elseif norm == "abs"
        return power
    else
        throw(ArgumentError("Unknown normalization method: $norm"))
    end
end

# Interface functions
freqs(ps::AbstractPowerSpectrum) = ps.freqs
power(ps::AbstractPowerSpectrum) = ps.power
errors(ps::AbstractPowerSpectrum) = ps.power_errors

# Base methods
Base.length(ps::AbstractPowerSpectrum) = length(ps.freqs)
Base.size(ps::AbstractPowerSpectrum) = (length(ps),)
Base.getindex(ps::AbstractPowerSpectrum, i) = (ps.freqs[i], ps.power[i], ps.power_errors[i])

# Show methods
function Base.show(io::IO, ps::PowerSpectrum)
    print(io, "PowerSpectrum(n=$(ps.n), df=$(ps.freqs[2]-ps.freqs[1]))")
end

function Base.show(io::IO, ps::AveragedPowerspectrum)
    print(io, "AveragedPowerspectrum(n=$(ps.n), m=$(ps.m), df=$(ps.freqs[2]-ps.freqs[1]))")
end

# Validation
function validate(ps::AbstractPowerSpectrum)
    if length(ps.freqs) != length(ps.power)
        throw(ArgumentError("Frequency and power arrays must have the same length"))
    end
    if length(ps.power) != length(ps.power_errors)
        throw(ArgumentError("Power and error arrays must have the same length"))
    end
    if any(x -> x < 0, ps.power)
        throw(ArgumentError("Power values must be non-negative"))
    end
    if any(x -> x < 0, ps.power_errors)
        throw(ArgumentError("Error values must be non-negative"))
    end
    return true
end