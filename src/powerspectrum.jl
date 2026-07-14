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

"""
    compute_rms(ps::Powerspectrum, min_freq::Real, max_freq::Real;
                poisson_noise_level=nothing)

Compute the fractional rms amplitude in the power spectrum between two
frequencies.

Mirrors Python Stingray's `Powerspectrum.compute_rms()`.

# Arguments
- `ps::Powerspectrum`: The power spectrum.
- `min_freq::Real`: Lower frequency bound.
- `max_freq::Real`: Upper frequency bound.

# Keyword Arguments
- `poisson_noise_level::Union{Nothing, Real}`: Poisson noise level in the same
  normalization as the power spectrum. If `nothing`, computed from the ideal
  Poisson level.

# Returns
- `(rms, rms_err)` — fractional rms amplitude and its uncertainty.
"""
function compute_rms(ps::Powerspectrum, min_freq::Real, max_freq::Real;
                     poisson_noise_level::Union{Nothing, Real}=nothing)

    good = (ps.freq .>= min_freq) .& (ps.freq .<= max_freq)

    # Handle per-bin k and m
    K_freq = ps.k isa AbstractVector ? ps.k[good] : ps.k
    M_freq = ps.m isa AbstractVector ? ps.m[good] : ps.m

    # Compute Poisson noise in unnormalized units
    if isnothing(poisson_noise_level)
        poisson_noise_unnorm = poisson_level("none"; n_ph=ps.nphots)
    else
        # Unnormalize the user-provided Poisson level
        poisson_noise_unnorm = _unnormalize_poisson(
            poisson_noise_level, ps.dt, ps.n, ps.nphots, ps.norm)
    end

    # Compute df for each frequency bin (accounting for rebinning)
    df_per_bin = if K_freq isa AbstractVector
        ps.df .* Float64.(K_freq)
    else
        ps.df * Float64(K_freq)
    end

    rms, rms_err = _get_rms_from_unnorm_periodogram(
        ps.unnorm_power[good], ps.nphots, df_per_bin;
        M=M_freq, poisson_noise_unnorm=poisson_noise_unnorm,
        segment_size=ps.segment_size
    )

    return rms, rms_err
end

"""
    _unnormalize_poisson(noise_level, dt, n, nphots, norm)

Convert a normalized Poisson noise level back to unnormalized units.
"""
function _unnormalize_poisson(noise_level::Real, dt::Real, n::Int,
                               nphots::Real, norm::String)
    mean_counts = nphots / n
    meanrate = mean_counts / dt

    if norm == "leahy"
        return noise_level * nphots / 2.0
    elseif norm == "frac"
        return noise_level * meanrate * nphots / 2.0
    elseif norm == "abs"
        return noise_level * nphots / (2.0 * meanrate)
    else  # "none"
        return noise_level
    end
end

"""
    _get_rms_from_unnorm_periodogram(unnorm_powers, nphots, df;
        M=1, poisson_noise_unnorm=nothing, segment_size=nothing)

Calculate fractional rms amplitude from unnormalized powers.

Mirrors Python Stingray's `fourier.get_rms_from_unnorm_periodogram()`.
"""
function _get_rms_from_unnorm_periodogram(
    unnorm_powers::AbstractVector, nphots::Real, df;
    M::Union{Int, Vector{Int}}=1,
    poisson_noise_unnorm::Union{Nothing, Real}=nothing,
    segment_size::Union{Nothing, Real}=nothing
)
    if isnothing(segment_size)
        segment_size = 1.0 / minimum(df isa AbstractVector ? df : [df])
    end

    if isnothing(poisson_noise_unnorm)
        poisson_noise_unnorm = nphots
    end

    meanrate = nphots / segment_size

    # Convert to fractional RMS normalization
    to_leahy(p) = p * 2.0 / nphots
    to_frac(p) = to_leahy(p) / meanrate

    powers_frac = to_frac.(unnorm_powers)
    poisson_frac = to_frac(poisson_noise_unnorm)

    # Compute total power error
    if M isa AbstractVector
        total_power_err = sqrt(sum(powers_frac .^ 2 .* (df isa AbstractVector ? df : fill(df, length(powers_frac))) .^ 2 ./ M))
    else
        df_arr = df isa AbstractVector ? df : fill(df, length(powers_frac))
        total_power_err = sqrt(sum(powers_frac .^ 2 .* df_arr .^ 2 ./ M))
    end

    # Subtract Poisson noise and integrate
    powers_sub = powers_frac .- poisson_frac
    df_arr = df isa AbstractVector ? df : fill(df, length(powers_frac))
    total_power_sub = sum(powers_sub .* df_arr)

    if total_power_sub < 0
        @warn "Poisson-subtracted power is below 0"
        return 0.0, 0.0
    end

    rms = sqrt(total_power_sub)
    rms_err = 0.5 / rms * total_power_err

    return rms, rms_err
end

