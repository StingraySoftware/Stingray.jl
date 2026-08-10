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

# ──────────────────────────────────────────────────────────────────────────────
# DynamicalPowerspectrum
# ──────────────────────────────────────────────────────────────────────────────

"""
    DynamicalPowerspectrum{T<:Real} <: AbstractPowerspectrum

A dynamical power spectrum (spectrogram) computed by dividing a single time
series into segments and computing the power spectrum for each.

The `dyn_ps` matrix has rows indexed by frequency and columns indexed by time
(segment midpoints). All values are real since it's an auto-spectrum.

Mirrors Python Stingray's `DynamicalPowerspectrum`.

# Fields
- `dyn_ps::Matrix{T}`: Frequency × time matrix of power spectral density.
- `freq::Vector{T}`: Mid-bin frequencies (Hz).
- `time::Vector{T}`: Mid-point time of each segment (s).
- `df::T`: Frequency resolution (Hz).
- `dt::T`: Time resolution of the spectrogram (= `segment_size` initially).
- `segment_size::T`: Length of each segment (s).
- `norm::String`: Normalization applied.
- `gti::Matrix{T}`: Good Time Intervals.
- `m::Int`: Number of averaged spectra per bin (1 initially).
- `nphots::T`: Mean photon count per segment.
- `meanrate::T`: Mean count rate (counts/s).
- `unnorm_conversion::T`: Mean ratio of normalized to unnormalized power.
- `type::String`: Always "dynamical_powerspectrum".
"""
struct DynamicalPowerspectrum{T<:Real} <: AbstractPowerspectrum
    dyn_ps::Matrix{T}
    freq::Vector{T}
    time::Vector{T}
    df::T
    dt::T
    segment_size::T
    norm::String
    gti::Matrix{T}
    m::Int
    nphots::T
    meanrate::T
    unnorm_conversion::T
    type::String
end

"""
    DynamicalPowerspectrum(ev::EventList, segment_size, dt; kwargs...)

Construct a `DynamicalPowerspectrum` from an `EventList`.

# Arguments
- `ev::EventList`: Event list containing photon arrival times.
- `segment_size::Real`: Length of each segment in seconds.
- `dt::Real`: Time resolution (sets Nyquist at 0.5/dt Hz).

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `silent::Bool=false`: Suppress progress bars.
"""
function DynamicalPowerspectrum(ev::EventList, segment_size::Real, dt::Real;
                                norm::String="frac",
                                silent::Bool=false)
    if segment_size < 2 * dt
        throw(ArgumentError("segment_size must be >= 2 * dt"))
    end

    ev_gti = has_gti(ev) ? gti(ev) : error("EventList must have GTIs")

    result = _pds_all_segments(
        times(ev), ev_gti,
        segment_size, dt;
        norm=norm, silent=silent
    )

    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    # Stack per-segment spectra into freq × time matrix
    dyn_ps = hcat(result.all_pds...)
    dyn_unnorm = hcat(result.all_unnorm_pds...)

    # Compute unnorm_conversion: mean(normalized / unnormalized)
    conv = dyn_ps ./ dyn_unnorm
    unnorm_conv = Float64(Statistics.mean(filter(!isnan, conv)))

    # Compute segment midpoints
    tstart, _ = time_intervals_from_gtis(ev_gti, segment_size)
    seg_times = Float64.(tstart .+ 0.5 * segment_size)

    n_bin = result.n_bin
    meanrate = result.nphots / n_bin / (segment_size / n_bin)

    return DynamicalPowerspectrum{Float64}(
        Float64.(dyn_ps),
        result.freq,
        seg_times,
        result.df,
        Float64(segment_size),
        Float64(segment_size),
        norm,
        ev_gti,
        1,
        result.nphots,
        meanrate,
        unnorm_conv,
        "dynamical_powerspectrum"
    )
end

"""
    DynamicalPowerspectrum(lc::LightCurve, segment_size; kwargs...)

Construct a `DynamicalPowerspectrum` from a `LightCurve`.

# Arguments
- `lc::LightCurve`: The input light curve.
- `segment_size::Real`: Length of each segment in seconds.

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `silent::Bool=false`: Suppress progress bars.
"""
function DynamicalPowerspectrum(lc::LightCurve, segment_size::Real;
                                norm::String="frac",
                                silent::Bool=false)
    dt_val = lc.dt isa AbstractVector ? lc.dt[1] : lc.dt

    if segment_size < 2 * dt_val
        throw(ArgumentError("segment_size must be >= 2 * dt"))
    end

    lc_gti = [lc.metadata.time_range[1] lc.metadata.time_range[2]]

    total_duration = lc_gti[1, 2] - lc_gti[1, 1]
    if segment_size > total_duration
        throw(ArgumentError(
            "segment_size ($segment_size s) is longer than total duration ($total_duration s)"))
    end

    result = _pds_all_segments(
        lc.time, lc_gti,
        segment_size, dt_val;
        norm=norm, silent=silent,
        fluxes=Float64.(lc.counts)
    )

    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    dyn_ps = hcat(result.all_pds...)
    dyn_unnorm = hcat(result.all_unnorm_pds...)
    conv = dyn_ps ./ dyn_unnorm
    unnorm_conv = Float64(Statistics.mean(filter(!isnan, conv)))

    tstart, _ = time_intervals_from_gtis(lc_gti, segment_size)
    seg_times = Float64.(tstart .+ 0.5 * segment_size)

    n_bin = result.n_bin
    meanrate = result.nphots / n_bin / dt_val

    return DynamicalPowerspectrum{Float64}(
        Float64.(dyn_ps),
        result.freq,
        seg_times,
        result.df,
        Float64(segment_size),
        Float64(segment_size),
        norm,
        lc_gti,
        1,
        result.nphots,
        meanrate,
        unnorm_conv,
        "dynamical_powerspectrum"
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# DynamicalPowerspectrum methods
# ──────────────────────────────────────────────────────────────────────────────

"""
    rebin_frequency(dps::DynamicalPowerspectrum, df_new; method=:mean)

Rebin the frequency axis of the dynamical power spectrum to a new resolution.

# Arguments
- `dps::DynamicalPowerspectrum`: The dynamical power spectrum.
- `df_new::Real`: New frequency resolution (must be >= current `df`).

# Keyword Arguments
- `method::Symbol=:mean`: Rebinning method (`:mean` or `:sum`).

# Returns
A new `DynamicalPowerspectrum` with the rebinned frequency axis.
"""
function rebin_frequency(dps::DynamicalPowerspectrum, df_new::Real; method::Symbol=:mean)
    n_time = size(dps.dyn_ps, 2)

    # Rebin each time column
    new_freq, col1, _, _ = rebin_data(dps.freq, dps.dyn_ps[:, 1], df_new; method=method)
    n_freq_new = length(col1)

    new_dyn = Matrix{Float64}(undef, n_freq_new, n_time)
    for j in 1:n_time
        _, ybin, _, _ = rebin_data(dps.freq, dps.dyn_ps[:, j], df_new; method=method)
        new_dyn[:, j] = ybin
    end

    return DynamicalPowerspectrum{Float64}(
        new_dyn, Float64.(new_freq), dps.time,
        Float64(df_new), dps.dt, dps.segment_size,
        dps.norm, dps.gti, dps.m,
        dps.nphots, dps.meanrate, dps.unnorm_conversion,
        dps.type
    )
end

"""
    rebin_time(dps::DynamicalPowerspectrum, dt_new; method=:mean)

Rebin the time axis of the dynamical power spectrum to a new resolution.

# Arguments
- `dps::DynamicalPowerspectrum`: The dynamical power spectrum.
- `dt_new::Real`: New time resolution (must be >= current `dt`).

# Keyword Arguments
- `method::Symbol=:mean`: Rebinning method (`:mean` or `:sum`).

# Returns
A new `DynamicalPowerspectrum` with the rebinned time axis.
"""
function rebin_time(dps::DynamicalPowerspectrum, dt_new::Real; method::Symbol=:mean)
    n_freq = size(dps.dyn_ps, 1)

    new_time, row1, _, _ = rebin_data(dps.time, dps.dyn_ps[1, :], dt_new; method=method)
    n_time_new = length(row1)

    new_dyn = Matrix{Float64}(undef, n_freq, n_time_new)
    for i in 1:n_freq
        _, ybin, _, _ = rebin_data(dps.time, dps.dyn_ps[i, :], dt_new; method=method)
        new_dyn[i, :] = ybin
    end

    return DynamicalPowerspectrum{Float64}(
        new_dyn, dps.freq, Float64.(new_time),
        dps.df, Float64(dt_new), dps.segment_size,
        dps.norm, dps.gti, dps.m,
        dps.nphots, dps.meanrate, dps.unnorm_conversion,
        dps.type
    )
end

"""
    rebin_by_n_intervals(dps::DynamicalPowerspectrum, n; method=:mean)

Average `n` consecutive time segments of the dynamical power spectrum.

# Arguments
- `dps::DynamicalPowerspectrum`: The dynamical power spectrum.
- `n::Int`: Number of consecutive segments to average.

# Keyword Arguments
- `method::Symbol=:mean`: `:mean` (average) or `:sum`.

# Returns
A new `DynamicalPowerspectrum` with fewer time columns and updated `m` and `dt`.
"""
function rebin_by_n_intervals(dps::DynamicalPowerspectrum, n::Int; method::Symbol=:mean)
    n_freq, n_time = size(dps.dyn_ps)
    n_new = n_time ÷ n

    if n_new < 1
        throw(ArgumentError("n=$n is larger than the number of time segments ($n_time)"))
    end

    new_dyn = Matrix{Float64}(undef, n_freq, n_new)
    new_time = Vector{Float64}(undef, n_new)

    for j in 1:n_new
        cols = ((j-1)*n + 1):(j*n)
        if method == :mean
            new_dyn[:, j] = Statistics.mean(dps.dyn_ps[:, cols], dims=2)
        else
            new_dyn[:, j] = sum(dps.dyn_ps[:, cols], dims=2)
        end
        new_time[j] = Statistics.mean(dps.time[cols])
    end

    new_m = dps.m * n

    return DynamicalPowerspectrum{Float64}(
        new_dyn, dps.freq, new_time,
        dps.df, dps.dt * n, dps.segment_size,
        dps.norm, dps.gti, new_m,
        dps.nphots, dps.meanrate, dps.unnorm_conversion,
        dps.type
    )
end

"""
    trace_maximum(dps::DynamicalPowerspectrum; min_freq=nothing, max_freq=nothing)

Find the frequency index of the maximum power for each time segment.

# Arguments
- `dps::DynamicalPowerspectrum`: The dynamical power spectrum.

# Keyword Arguments
- `min_freq::Union{Nothing, Real}=nothing`: Lower frequency bound.
- `max_freq::Union{Nothing, Real}=nothing`: Upper frequency bound.

# Returns
A `Vector{Int}` of frequency indices (into `dps.freq`) of the peak power
for each time segment.
"""
function trace_maximum(dps::DynamicalPowerspectrum;
                       min_freq::Union{Nothing, Real}=nothing,
                       max_freq::Union{Nothing, Real}=nothing)
    freq = dps.freq
    mask = trues(length(freq))
    if !isnothing(min_freq)
        mask .&= freq .>= min_freq
    end
    if !isnothing(max_freq)
        mask .&= freq .<= max_freq
    end
    idx_range = findall(mask)

    n_time = size(dps.dyn_ps, 2)
    result = Vector{Int}(undef, n_time)
    for j in 1:n_time
        col = dps.dyn_ps[idx_range, j]
        local_idx = argmax(col)
        result[j] = idx_range[local_idx]
    end
    return result
end
