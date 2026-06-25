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

    spectral = _extract_spectral_result(result; complex_power=false)
    params = _compute_spectral_params(segment_size, dt, ev_gti)
    power_err, unnorm_power_err = _compute_power_errors(
        spectral.power, spectral.unnorm_power, params.m)
    stats = _compute_photon_stats(length(times(ev)), params.m, segment_size)

    return Powerspectrum{Float64}(
        spectral.freq, spectral.power, power_err,
        spectral.unnorm_power, unnorm_power_err,
        params.df, params.dt, params.n_bin, params.m,
        1,  # k = 1 (not rebinned)
        stats.nphots, norm,
        ev_gti, Float64(segment_size),
        nothing, "poisson", "powerspectrum"
    )
end

"""
    Powerspectrum(lc::LightCurve, segment_size::Real; kwargs...)

Construct a `Powerspectrum` from a `LightCurve`.

Extracts times, counts, and derives GTI from the light curve's time range.
The dt is taken from `lc.dt` (the bin size stored on the LightCurve struct).

# Arguments
- `lc::LightCurve`: The input light curve.
- `segment_size::Real`: Length of each segment in seconds.

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `use_common_mean::Bool=true`: Use the mean from the full light curve.
- `silent::Bool=false`: Suppress progress bars.
"""
function Powerspectrum(lc::LightCurve, segment_size::Real;
                       norm::String="frac",
                       use_common_mean::Bool=true,
                       silent::Bool=false)

    # dt is a direct field on LightCurve, handle Union{T,Vector{T}}
    dt_val = lc.dt isa AbstractVector ? lc.dt[1] : lc.dt

    # Derive GTI from light curve time range (LightCurveMetadata has no gti field)
    lc_gti = [lc.metadata.time_range[1] lc.metadata.time_range[2]]

    result = avg_pds_from_events(
        lc.time, lc_gti,
        segment_size, dt_val;
        norm=norm,
        use_common_mean=use_common_mean,
        silent=silent,
        fluxes=Float64.(lc.counts)
    )

    spectral = _extract_spectral_result(result; complex_power=false)
    params = _compute_spectral_params(segment_size, dt_val, lc_gti)
    power_err, unnorm_power_err = _compute_power_errors(
        spectral.power, spectral.unnorm_power, params.m)
    stats = _compute_photon_stats(sum(lc.counts), params.m, segment_size)

    err_dist = lc.err_method == :gaussian ? "gauss" : "poisson"
    variance_val = err_dist == "gauss" ? Float64(var(lc.counts)) : nothing

    return Powerspectrum{Float64}(
        spectral.freq, spectral.power, power_err,
        spectral.unnorm_power, unnorm_power_err,
        params.df, params.dt, params.n_bin, params.m,
        1,
        stats.nphots, norm,
        lc_gti, Float64(segment_size),
        variance_val, err_dist, "powerspectrum"
    )
end
