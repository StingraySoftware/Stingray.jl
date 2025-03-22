using DSP
using Statistics
using LinearAlgebra

# Base Powerspectrum type
struct Powerspectrum
    freq::Vector{Float64}
    power::Vector{Float64}
    unnorm_power::Vector{Float64}
    df::Float64
    dt::Float64
    nphots::Int
    m::Int
    n::Int
    k::Int
    norm::String
end

# Constructor for Powerspectrum
function Powerspectrum(data::Vector{Float64}, dt::Float64, norm::String="frac")
    n = length(data)
    freq = rfftfreq(n, 1.0 / dt)
    power = abs2.(rfft(data)) / n
    unnorm_power = copy(power)
    
    # Normalization
    if norm == "leahy"
        power *= 2.0 / sum(data)
    elseif norm == "frac"
        power /= sum(data)^2
    end

    Powerspectrum(freq, power, unnorm_power, freq[2] - freq[1], dt, sum(data), 1, n, 1, norm)
end

# Rebin the power spectrum
function rebin(ps::Powerspectrum, df::Float64)
    new_freq = ps.freq[1]:df:ps.freq[end]
    new_power = zeros(length(new_freq))
    for i in 1:length(new_freq)
        idx = (ps.freq .>= new_freq[i]) .& (ps.freq .< new_freq[i] + df)
        new_power[i] = mean(ps.power[idx])
    end
    Powerspectrum(new_freq, new_power, ps.unnorm_power, df, ps.dt, ps.nphots, ps.m, ps.n, ps.k, ps.norm)
end

# Compute RMS in a frequency range
function compute_rms(ps::Powerspectrum, min_freq::Float64, max_freq::Float64)
    idx = (ps.freq .>= min_freq) .& (ps.freq .<= max_freq)
    rms = sqrt(sum(ps.power[idx]) * ps.df)
    rms
end

# Averaged Powerspectrum type
struct AveragedPowerspectrum
    freq::Vector{Float64}
    power::Vector{Float64}
    unnorm_power::Vector{Float64}
    df::Float64
    dt::Float64
    nphots::Int
    m::Int
    n::Int
    k::Int
    norm::String
end

# Constructor for AveragedPowerspectrum
function AveragedPowerspectrum(data::Vector{Vector{Float64}}, dt::Float64, norm::String="frac")
    n_segments = length(data)
    n = length(data[1])
    freq = rfftfreq(n, 1.0 / dt)
    power = zeros(length(freq))
    unnorm_power = zeros(length(freq))
    
    for segment in data
        ft = rfft(segment)
        power .+= abs2.(ft) / n
        unnorm_power .+= abs2.(ft) / n
    end
    power ./= n_segments
    unnorm_power ./= n_segments

    # Normalization
    if norm == "leahy"
        power .*= 2.0 / mean(sum.(data))
    elseif norm == "frac"
        power ./= mean(sum.(data))^2
    end

    AveragedPowerspectrum(freq, power, unnorm_power, freq[2] - freq[1], dt, sum(sum.(data)), n_segments, n, 1, norm)
end

# Dynamical Powerspectrum type
struct DynamicalPowerspectrum
    time::Vector{Float64}
    freq::Vector{Float64}
    dyn_ps::Matrix{Float64}
    df::Float64
    dt::Float64
    nphots::Int
    m::Int
    norm::String
end

# Constructor for DynamicalPowerspectrum
function DynamicalPowerspectrum(lc::Vector{Float64}, dt::Float64, segment_size::Float64, norm::String="frac")
    n_segments = Int(floor(length(lc) * dt / segment_size))
    n = Int(segment_size / dt)
    freq = rfftfreq(n, 1.0 / dt)
    dyn_ps = zeros(length(freq), n_segments)
    time = zeros(n_segments)
    
    for i in 1:n_segments
        start_idx = Int((i - 1) * n + 1)
        end_idx = Int(i * n)
        segment = lc[start_idx:end_idx]
        ft = rfft(segment)
        dyn_ps[:, i] = abs2.(ft) / n
        time[i] = (start_idx + end_idx) / 2 * dt
    end

    # Normalization
    if norm == "leahy"
        dyn_ps .*= 2.0 / mean(sum.(eachcol(dyn_ps)))
    elseif norm == "frac"
        dyn_ps ./= mean(sum.(eachcol(dyn_ps)))^2
    end

    DynamicalPowerspectrum(time, freq, dyn_ps, freq[2] - freq[1], segment_size, sum(lc), 1, norm)
end

# Example usage
function main()
    # Create a synthetic light curve
    dt = 0.1
    times = 0:dt:100
    counts = 100 .+ 10 * sin.(2 * π * 0.5 * times)

    # Compute a single power spectrum
    ps = Powerspectrum(counts, dt)
    println("Single Power Spectrum Frequencies: ", ps.freq)
    println("Single Power Spectrum Powers: ", ps.power)

    # Compute an averaged power spectrum
    segment_size = 10.0
    n_segments = Int(floor(length(counts) * dt / segment_size))
    segments = [counts[Int((i-1)*segment_size/dt+1):Int(i*segment_size/dt)] for i in 1:n_segments]
    aps = AveragedPowerspectrum(segments, dt)
    println("Averaged Power Spectrum Frequencies: ", aps.freq)
    println("Averaged Power Spectrum Powers: ", aps.power)

    # Compute a dynamical power spectrum
    dps = DynamicalPowerspectrum(counts, dt, segment_size)
    println("Dynamical Power Spectrum Time: ", dps.time)
    println("Dynamical Power Spectrum Frequencies: ", dps.freq)
    println("Dynamical Power Spectrum Matrix: ", dps.dyn_ps)
end

main()