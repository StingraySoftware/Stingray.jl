abstract type AbstractPowerSpectrum end

struct PowerSpectrum <: AbstractPowerSpectrum
    freqs::Vector{Float64}
    power::Vector{Float64}
    power_errors::Vector{Float64}
    dt::Float64
    n::Int
end

struct AveragedPowerspectrum <: AbstractPowerSpectrum
    freqs::Vector{Float64}
    power::Vector{Float64}
    power_errors::Vector{Float64}
    m::Int
    n::Int
    dt::Float64
end

# Power spectrum functions
function powerspectrum(events::EventList{T}; dt::Real = 0.001,
                      segment_size::Int = length(events.times) ÷ 32, 
                      norm::String = "frac") where T
    @assert segment_size > 0 "Segment size must be positive"
    @assert length(events.times) >= segment_size "Not enough data points"
    return _powerspectrum_averaged(events, float(dt), segment_size, norm)
end

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

# Helper functions
function hanning(N::Int)
    return [0.5 * (1 - cos(2π * n / (N-1))) for n in 0:N-1]
end

function _prepare_segments(signal::Vector{Float64}, segment_size::Int)
    n_segments = length(signal) ÷ segment_size
    segments = [signal[(i-1)*segment_size+1 : i*segment_size] for i in 1:n_segments]
    return segments
end

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