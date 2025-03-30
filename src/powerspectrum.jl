"""
Abstract type representing a power spectrum, which characterizes the distribution 
of power across different frequencies in a signal.

Subtypes include:
- PowerSpectrum{T}: Represents a power spectrum for a single signal segment
- AveragedPowerspectrum{T}: Represents a power spectrum averaged over multiple segments
"""
abstract type AbstractPowerSpectrum{T} end

"""
Represents a power spectrum computed from a single signal segment.

## Fields
- `freqs::Vector{T}`: Vector of frequency values 
- `power::Vector{T}`: Corresponding power values for each frequency
- `power_errors::Vector{T}`: Uncertainty estimates for power at each frequency
- `dt::T`: Time sampling interval
- `n::Int`: Number of data points in the original segment
"""
struct PowerSpectrum{T} <: AbstractPowerSpectrum{T}
    freqs::Vector{T}
    power::Vector{T}
    power_errors::Vector{T}
    dt::T
    n::Int
end

"""
Represents a power spectrum computed by averaging multiple signal segments.

## Fields
- `freqs::Vector{T}`: Vector of frequency values
- `power::Vector{T}`: Average power values for each frequency
- `power_errors::Vector{T}`: Standard error of power estimates
- `m::Int`: Number of segments used in averaging
- `n::Int`: Number of data points in each segment
- `dt::T`: Time sampling interval
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
    powerspectrum(events::EventList{T}; dt::Real=0.001,
                 segment_size::Int=length(events.times) ÷ 32, 
                 norm::Symbol=:frac) where T -> AveragedPowerspectrum{T}

Compute the power spectrum for a list of events using segmented averaging.

## Arguments
- `events::EventList{T}`: Input event list containing times and optional energies
- `dt::Real=0.001`: Time sampling interval (default: 0.001)
- `segment_size::Int=length(events.times) ÷ 32`: Size of each segment for FFT (default: 1/32 of total data)
- `norm::Symbol=:frac`: Normalization method for power spectrum
    - `:leahy`: Multiplies power by 2
    - `:frac`: Normalizes power by mean power
    - `:abs`: No normalization

## Returns
- `AveragedPowerspectrum{T}`: Power spectrum computed from segmented data

## Throws
- `ArgumentError` if segment size is invalid or insufficient data points

## Examples
```julia
events = EventList(times, energies, metadata)
ps = powerspectrum(events, dt=0.01, segment_size=1024, norm=:leahy)
```
"""
function powerspectrum(events::EventList{T}; dt::Real = 0.001,
                      segment_size::Int = length(events.times) ÷ 32, 
                      norm::Symbol = :frac) where T
    if segment_size <= 0
        throw(ArgumentError("Segment size must be positive"))
    end
    if length(events.times) < segment_size
        throw(ArgumentError("Not enough data points"))
    end
    return _powerspectrum_averaged(events, convert(T, dt), segment_size, norm)
end

"""
    _powerspectrum_averaged(events::EventList{T}, dt::T, 
                          segment_size::Int, norm::Symbol) where T -> AveragedPowerspectrum{T}

Compute averaged power spectrum from event list segments.

## Arguments
- `events::EventList{T}`: Input event list
- `dt::T`: Time sampling interval
- `segment_size::Int`: Number of points in each segment
- `norm::Symbol`: Normalization method for power spectrum

## Returns
- `AveragedPowerspectrum{T}`: Power spectrum averaged over multiple segments

## Notes
- Prioritizes using event energies if available, otherwise uses event times
- Computes FFT for each segment and averages power
- Calculates frequency-dependent power with specified normalization
"""
function _powerspectrum_averaged(events::EventList{T}, 
                               dt::T, 
                               segment_size::Int, 
                               norm::Symbol) where T
    signal = events.times
    
    if !all(iszero, events.energies)
        signal = events.energies
    end
    
    n_segments = length(signal) ÷ segment_size
    segments = [signal[(i-1)*segment_size+1 : i*segment_size] for i in 1:n_segments]
    
    # Apply Hanning window if needed
    windowed_segments = [segment .* hanning(segment_size) for segment in segments]
    
    ffts = [fft(segment) for segment in windowed_segments]
    
    nyquist_index = segment_size ÷ 2 + 1
    powers = [abs2.(fft[1:nyquist_index]) for fft in ffts]
    
    avg_power = mean(powers)
    
    freqs = FFTW.rfftfreq(segment_size, 1/dt)
    
    if length(freqs) != length(avg_power)
        throw(ArgumentError("Frequency and power arrays must match"))
    end
    
    normalized_power = _normalize_power(avg_power, norm)
    power_errors = std(powers) ./ sqrt(length(segments))
    
    positive_indices = 2:length(freqs)
    
    return AveragedPowerspectrum{T}(
        convert(Vector{T}, freqs[positive_indices]),
        convert(Vector{T}, normalized_power[positive_indices]),
        convert(Vector{T}, power_errors[positive_indices]),
        length(segments),
        segment_size,
        dt
    )
end

"""
    hanning(N::Int) -> Vector{Float64}

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
    _normalize_power(power::AbstractVector, norm::Symbol) -> Vector{Float64}

Normalize power spectrum based on specified method.

## Arguments
- `power`: Input power values
- `norm`: Normalization method (`:leahy`, `:frac`, or `:abs`)

## Returns
- `Vector{Float64}`: Normalized power values

## Throws
- `ArgumentError` for unknown normalization method
"""
function _normalize_power(power::AbstractVector, norm::Symbol)
    if norm == :leahy
        return power .* 2
    elseif norm == :frac
        return power ./ mean(power)
    elseif norm == :abs
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
function Base.show(io::IO, ps::PowerSpectrum{T}) where T
    print(io, "PowerSpectrum{$T}(n=$(ps.n), df=$(ps.freqs[2]-ps.freqs[1]))")
end

function Base.show(io::IO, ps::AveragedPowerspectrum{T}) where T
    print(io, "AveragedPowerspectrum{$T}(n=$(ps.n), m=$(ps.m), df=$(ps.freqs[2]-ps.freqs[1]))")
end

"""
    validate(ps::AbstractPowerSpectrum) -> Bool

Validate the power spectrum structure.

## Returns
- `true` if valid, throws ArgumentError otherwise

## Throws
- `ArgumentError` if any validation checks fail
"""
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