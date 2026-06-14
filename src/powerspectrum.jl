"""
    Powerspectrum{T<:Real} <: AbstractPowerspectrum

A power spectrum (periodogram) computed from a single time series.

A power spectrum is mathematically a cross-spectrum of a signal with itself,
which is why `AbstractPowerspectrum <: AbstractCrossspectrum`. Unlike the
`Crossspectrum`, the power values are real-valued (the imaginary part of a
signal crossed with itself is always zero).

# Fields
- `freq::Vector{T}`: Mid-bin frequencies of the Fourier transform.
- `power::Vector{T}`: Normalized power spectral density (real-valued).
- `power_err::Vector{T}`: Uncertainties on `power`, approximated as `power / √m`.
- `unnorm_power::Vector{T}`: Unnormalized powers.
- `unnorm_power_err::Vector{T}`: Uncertainties on `unnorm_power`.
- `df::T`: Frequency resolution (Hz).
- `dt::T`: Time resolution of the input data (s).
- `n::Int`: Number of bins per segment.
- `m::Int`: Number of averaged power spectra.
- `k::Union{Int, Vector{Int}}`: Rebinning factor (1 if not rebinned).
- `nphots::T`: Total number of photons in the light curve.
- `norm::String`: Normalization applied: "frac", "abs", "leahy", or "none".
- `gti::Matrix{T}`: Good Time Intervals, Nx2 matrix of [start, stop] pairs.
- `segment_size::T`: Length of each segment in seconds.
- `variance::Union{Nothing, T}`: Variance of the light curve (if Gaussian errors were used).
- `err_dist::String`: Error distribution assumption: "poisson" or "gauss".
- `type::String`: Object type identifier, always "powerspectrum".
"""
struct Powerspectrum{T<:Real} <: AbstractPowerspectrum
    freq::Vector{T}
    power::Vector{T}
    power_err::Vector{T}
    unnorm_power::Vector{T}
    unnorm_power_err::Vector{T}
    df::T
    dt::T
    n::Int
    m::Int
    k::Union{Int, Vector{Int}}
    nphots::T
    norm::String
    gti::Matrix{T}
    segment_size::T
    variance::Union{Nothing, T}
    err_dist::String
    type::String
end

function Base.show(io::IO, ps::Powerspectrum{T}) where T
    print(io, "Powerspectrum($(ps.norm), $(ps.m) segments, ",
          "$(length(ps.freq)) freq bins, ",
          "df=$(round(ps.df, sigdigits=4)) Hz, ",
          "freq range=[$(round(ps.freq[1], sigdigits=4)), ",
          "$(round(ps.freq[end], sigdigits=4))] Hz)")
end

"""
    Powerspectrum(ev::EventList, segment_size::Real, dt::Real; kwargs...)

Construct a `Powerspectrum` from an `EventList`.

Calls `avg_pds_from_events` and computes metadata from input parameters.
Error bars are computed as `power / √m`.

# Arguments
- `ev::EventList`: The event list containing photon arrival times.
- `segment_size::Real`: Length of each segment in seconds.
- `dt::Real`: Time resolution (sets Nyquist frequency at 0.5/dt Hz).

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `use_common_mean::Bool=true`: Use the mean from the full light curve.
- `silent::Bool=false`: Suppress progress bars.
"""
function Powerspectrum(ev::EventList, segment_size::Real, dt::Real;
                       norm::String="frac",
                       use_common_mean::Bool=true,
                       silent::Bool=false)

    ev_gti = has_gti(ev) ? gti(ev) : error("EventList must have GTIs")

    result = avg_pds_from_events(
        times(ev), ev_gti,
        segment_size, dt;
        norm=norm,
        use_common_mean=use_common_mean,
        silent=silent
    )

    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    # Extract spectral data from the DataFrame
    freq = Float64.(result[!, "freq"])
    power = Float64.(real.(result[!, "power"]))
    unnorm_power = Float64.(real.(result[!, "unnorm_power"]))

    # Compute metadata from input parameters
    n_bin = round(Int, segment_size / dt)
    dt_adj = segment_size / n_bin
    df_val = 1.0 / segment_size
    m_val = _count_valid_segments(ev_gti, segment_size)

    power_err = power ./ sqrt(m_val)
    unnorm_power_err = unnorm_power ./ sqrt(m_val)

    nphots_val = Float64(length(times(ev))) / m_val

    return Powerspectrum{Float64}(
        freq, power, power_err, unnorm_power, unnorm_power_err,
        df_val, dt_adj, n_bin, m_val,
        1,  # k = 1 (not rebinned)
        nphots_val, norm,
        ev_gti, Float64(segment_size),
        nothing, "poisson", "powerspectrum"
    )
end
