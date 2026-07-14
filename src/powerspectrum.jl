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
- `m::Union{Int, Vector{Int}}`: Number of averaged power spectra.
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
    m_display = ps.m isa Int ? ps.m : "$(length(ps.m))-element"
    print(io, "Powerspectrum($(ps.norm), $(m_display) segments, ",
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

# ──────────────────────────────────────────────────────────────────────────────
# Powerspectrum rebinning methods
# ──────────────────────────────────────────────────────────────────────────────


"""
    rebin(ps::Powerspectrum{T}, df_new::Real; f=nothing, method=:mean)

Rebin the power spectrum to a new frequency resolution `df_new`.

Mirrors Python Stingray's `Powerspectrum.rebin()`, which delegates to
`Crossspectrum.rebin()` and sets `nphots = nphots1`.

# Arguments
- `ps::Powerspectrum{T}`: The power spectrum to rebin.
- `df_new::Real`: The new frequency resolution.

# Keyword Arguments
- `f::Union{Nothing, Real}`: Rebin factor. If specified, `df_new = f * ps.df`.
- `method::Symbol`: `:mean` or `:sum`. Default `:mean`.

# Returns
A new `Powerspectrum` with updated `freq`, `power`, `df`, `m`, and `k`.
"""
function rebin(ps::Powerspectrum{T}, df_new::Real;
               f::Union{Nothing, Real}=nothing,
               method::Symbol=:mean) where T

    if !isnothing(f)
        df_new = f * ps.df
    end

    # Rebin normalized power
    binfreq, binpower, binerr, step_size = rebin_data(
        ps.freq, ps.power, df_new;
        yerr=ps.power_err, method=method, dx=ps.df
    )

    # Rebin unnormalized power
    _, binpower_unnorm, binpower_err_unnorm, _ = rebin_data(
        ps.freq, ps.unnorm_power, df_new;
        yerr=ps.unnorm_power_err, method=method, dx=ps.df
    )

    # Update m: m_new = round(Int, step_size * m_old)
    m_old = ps.m isa Int ? ps.m : round(Int, mean(ps.m))
    new_m = [round(Int, s * m_old) for s in step_size]
    new_m_val = length(unique(new_m)) == 1 ? new_m[1] : new_m

    return Powerspectrum{Float64}(
        Float64.(binfreq),
        Float64.(real.(binpower)),
        Float64.(real.(binerr)),
        Float64.(real.(binpower_unnorm)),
        Float64.(real.(binpower_err_unnorm)),
        Float64(df_new), ps.dt, ps.n, new_m_val,
        round(Int, mean(step_size)),  # k = rebinning factor
        ps.nphots, ps.norm,
        ps.gti, ps.segment_size,
        ps.variance, ps.err_dist, ps.type
    )
end

"""
    rebin_log(ps::Powerspectrum{T}, f::Real=0.01)

Logarithmically rebin the power spectrum. Each new frequency bin grows as
`dν_j = dν_{j-1} * (1+f)`.

Mirrors Python Stingray's `Powerspectrum.rebin_log()` via `Crossspectrum.rebin_log()`.

# Arguments
- `ps::Powerspectrum{T}`: The power spectrum to rebin.
- `f::Real`: Growth factor for frequency bins. Default `0.01`.

# Returns
A new `Powerspectrum` with updated `freq`, `power`, `m` (vector), `k` (vector).
"""
function rebin_log(ps::Powerspectrum{T}, f::Real=0.01) where T

    # Rebin normalized power
    binfreq, binpower, binpower_err, nsamples = rebin_data_log(
        ps.freq, ps.power, f;
        y_err=ps.power_err, dx=ps.df
    )

    # Rebin unnormalized power
    _, binpower_unnorm, binpower_err_unnorm, _ = rebin_data_log(
        ps.freq, ps.unnorm_power, f;
        y_err=ps.unnorm_power_err, dx=ps.df
    )

    # m = nsamples * original_m
    m_old = ps.m isa Int ? ps.m : round(Int, mean(ps.m))
    new_m = nsamples .* m_old

    # df = median of new frequency spacing
    new_df = length(binfreq) > 1 ? Float64(median(diff(binfreq))) : ps.df

    return Powerspectrum{Float64}(
        Float64.(binfreq),
        Float64.(real.(binpower)),
        Float64.(real.(binpower_err)),
        Float64.(real.(binpower_unnorm)),
        Float64.(real.(binpower_err_unnorm)),
        new_df, ps.dt, ps.n, new_m,
        nsamples,  # k = nsamples (Vector{Int})
        ps.nphots, ps.norm,
        ps.gti, ps.segment_size,
        ps.variance, ps.err_dist, ps.type
    )
end
