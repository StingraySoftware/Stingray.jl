"""
Abstract base type for power spectrum implementations.
Created: 2025-03-25 23:59:59 UTC
Author: kashish2210
"""
abstract type AbstractPowerSpectrum end

"""
    PowerSpectrum

Basic power spectrum type.

## Fields
- `freqs::Vector{Float64}`: Frequency bins
- `power::Vector{Float64}`: Power values
- `power_errors::Vector{Float64}`: Uncertainties in power values
- `dt::Float64`: Time resolution
- `n::Int`: Number of points
"""
struct PowerSpectrum <: AbstractPowerSpectrum
    freqs::Vector{Float64}
    power::Vector{Float64}
    power_errors::Vector{Float64}
    dt::Float64
    n::Int
end

"""
    AveragedPowerspectrum <: AbstractPowerSpectrum

Power spectrum averaged over multiple segments.

## Fields
- `freqs::Vector{Float64}`: Frequency bins
- `power::Vector{Float64}`: Power values
- `power_errors::Vector{Float64}`: Uncertainties in power values
- `m::Int`: Number of segments averaged
- `n::Int`: Number of points per segment
- `dt::Float64`: Time resolution
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
    powerspectrum(events::AbstractEventList; kwargs...)

Compute power spectrum from event data.

## Parameters
- `dt::Real`: Time resolution
- `segment_size::Int=0`: Size of segments (0 for no averaging)
- `norm::String="frac"`: Normalization ("frac", "leahy", or "abs")

## Returns
- `PowerSpectrum` or `AveragedPowerspectrum`
"""
function powerspectrum(events::AbstractEventList{T};
                      dt::Real,
                      segment_size::Int=0,
                      norm::String="frac") where T
    
    if segment_size == 0
        return _powerspectrum_single(events, dt, norm)
    else
        return _powerspectrum_averaged(events, dt, segment_size, norm)
    end
end

# Internal implementation functions
function _powerspectrum_single(events::AbstractEventList{T}, dt::Real, norm::String) where T
    # Validation
    if dt <= 0
        throw(ArgumentError("dt must be positive"))
    end
    
    event_times = times(events)
    
    # Create time bins
    tstart = minimum(event_times)
    tend = maximum(event_times)
    edges = range(tstart, tend, step=dt)
    
    # Bin events
    counts = fit(Histogram, event_times, edges).weights
    
    # Compute FFT
    ft = rfft(counts)
    freqs = FFTW.rfftfreq(length(counts), 1/dt)
    power = abs2.(ft)
    
    # Apply normalization
    power, errors = _normalize_power(power, counts, norm)
    
    PowerSpectrum(
        freqs[2:end],
        power[2:end],
        errors[2:end],
        dt,
        length(counts)
    )
end

function _powerspectrum_averaged(events::AbstractEventList{T}, 
                               dt::Real, 
                               segment_size::Int, 
                               norm::String) where T
    # Validation
    if dt <= 0
        throw(ArgumentError("dt must be positive"))
    end
    if segment_size <= 0
        throw(ArgumentError("segment_size must be positive"))
    end
    
    event_times = times(events)
    
    # Calculate segments
    tstart = minimum(event_times)
    tend = maximum(event_times)
    total_time = tend - tstart
    n_segments = floor(Int, total_time / (segment_size * dt))
    
    if n_segments == 0
        throw(ArgumentError("Time series too short for given segment_size"))
    end
    
    # Initialize arrays
    freqs = FFTW.rfftfreq(segment_size, 1/dt)
    powers = zeros(Float64, length(freqs), n_segments)
    
    # Process segments
    for i in 1:n_segments
        t_start = tstart + (i-1) * segment_size * dt
        t_end = t_start + segment_size * dt
        
        # Select events in this time window
        mask = (event_times .>= t_start) .& (event_times .< t_end)
        segment_times = event_times[mask]
        segment_energies = energies(events)[mask]
        
        # Create segment event list
        segment_events = EventList(
            events.filename,
            segment_times,
            segment_energies,
            events.metadata
        )
        
        segment_spectrum = _powerspectrum_single(segment_events, dt, norm)
        powers[2:end, i] = segment_spectrum.power
    end
    
    # Average results
    avg_power = mean(powers, dims=2)[:, 1]
    power_err = std(powers, dims=2)[:, 1] / sqrt(n_segments)
    
    AveragedPowerspectrum(
        freqs[2:end],
        avg_power[2:end],
        power_err[2:end],
        n_segments,
        segment_size,
        dt
    )
end

function _normalize_power(power::Vector{Float64}, 
                        counts::Vector{Float64}, 
                        norm::String)
    if norm == "leahy"
        mean_rate = mean(counts)
        normalized = 2 .* power ./ mean_rate
        errors = 2 .* sqrt.(power) ./ mean_rate
    elseif norm == "frac"
        mean_rate = mean(counts)
        normalized = power ./ (mean_rate^2)
        errors = sqrt.(power) ./ mean_rate^2
    elseif norm == "abs"
        normalized = power
        errors = sqrt.(power)
    else
        throw(ArgumentError("Unknown normalization: $norm"))
    end
    
    return normalized, errors
end

"""
    freqs(ps::AbstractPowerSpectrum)

Get the frequency bins of a power spectrum.
"""
freqs(ps::AbstractPowerSpectrum) = ps.freqs

"""
    power(ps::AbstractPowerSpectrum)

Get the power values of a power spectrum.
"""
power(ps::AbstractPowerSpectrum) = ps.power

"""
    errors(ps::AbstractPowerSpectrum)

Get the error estimates of a power spectrum.
"""
errors(ps::AbstractPowerSpectrum) = ps.power_errors

# Basic array interface methods
Base.length(ps::AbstractPowerSpectrum) = length(ps.freqs)
Base.size(ps::AbstractPowerSpectrum) = (length(ps),)
Base.getindex(ps::AbstractPowerSpectrum, i) = (ps.freqs[i], ps.power[i], ps.power_errors[i])

# Pretty printing
function Base.show(io::IO, ps::PowerSpectrum)
    print(io, "PowerSpectrum(n=$(ps.n), df=$(ps.freqs[2]-ps.freqs[1]))")
end

function Base.show(io::IO, ps::AveragedPowerspectrum)
    print(io, "AveragedPowerspectrum(n=$(ps.n), m=$(ps.m), df=$(ps.freqs[2]-ps.freqs[1]))")
end

# Validation function
"""
    validate(ps::AbstractPowerSpectrum)

Validate the power spectrum structure.
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