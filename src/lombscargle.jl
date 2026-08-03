# ──────────────────────────────────────────────────────────────────────────────
# Lomb-Scargle Periodograms for unevenly sampled data
#
# Mirrors Python Stingray's stingray/lombscargle.py
# ──────────────────────────────────────────────────────────────────────────────

"""
    autofrequency(; min_freq=nothing, max_freq=nothing, df=nothing,
                    dt=nothing, length=nothing, nyquist_factor=1)

Build the frequency grid for a Lomb-Scargle periodogram.

Determines `df`, `max_freq`, and `min_freq` from whatever subset of
parameters the user provides, then returns a regularly spaced frequency
array from `min_freq` to `max_freq` in steps of `df`.

Mirrors Python Stingray's `lombscargle._autofrequency`.

# Keyword Arguments
- `min_freq::Union{Nothing, Real}=nothing`: Minimum frequency. Defaults to `df / 2`.
- `max_freq::Union{Nothing, Real}=nothing`: Maximum frequency. Defaults to Nyquist: `nyquist_factor * 0.5 / dt`.
- `df::Union{Nothing, Real}=nothing`: Frequency resolution. Defaults to `1 / length`.
- `dt::Union{Nothing, Real}=nothing`: Time resolution of the data.
- `length::Union{Nothing, Real}=nothing`: Total duration of the observation.
- `nyquist_factor::Real=1`: Multiplier on the Nyquist frequency.

# Returns
- `Vector{Float64}`: The array of frequencies.

# Examples
```julia
freqs = autofrequency(min_freq=0.1, max_freq=0.5, df=0.1)
# [0.1, 0.2, 0.3, 0.4, 0.5]

freqs = autofrequency(min_freq=0.1, dt=1.0, length=10.0)
# [0.1, 0.2, 0.3, 0.4, 0.5]
```
"""
function autofrequency(; min_freq::Union{Nothing, Real}=nothing,
                         max_freq::Union{Nothing, Real}=nothing,
                         df::Union{Nothing, Real}=nothing,
                         dt::Union{Nothing, Real}=nothing,
                         length::Union{Nothing, Real}=nothing,
                         nyquist_factor::Real=1)

    # Determine df
    if (isnothing(df) || df <= 0) && isnothing(length)
        throw(ValueError("Either df or length must be specified."))
    elseif isnothing(df) || df <= 0
        df = 1.0 / length
    end

    # Determine max_freq
    if isnothing(max_freq) && (isnothing(dt) || dt == 0)
        throw(ValueError("Either max_freq or dt must be specified."))
    elseif isnothing(max_freq)
        max_freq = nyquist_factor * 0.5 / dt
    end

    # Determine min_freq
    if isnothing(min_freq)
        min_freq = df / 2.0
    elseif min_freq <= 0
        @warn "min_freq must be positive and >0. Setting to df / 2."
        min_freq = df / 2.0
    end

    freq = collect(Float64(min_freq):Float64(df):Float64(max_freq + df))
    # Trim any values beyond max_freq (floating point overshoot)
    filter!(f -> f <= max_freq + df * 0.5, freq)

    return freq
end

# ValueError alias for consistency with Python-style error messages
struct ValueError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ValueError) = print(io, "ValueError: ", e.msg)

# ──────────────────────────────────────────────────────────────────────────────
# LombScargleCrossspectrum
# ──────────────────────────────────────────────────────────────────────────────

"""
    LombScargleCrossspectrum{T<:Real} <: AbstractCrossspectrum

A cross spectrum computed from two unevenly sampled time series using
the Lomb-Scargle Fourier Transform.

Unlike the standard `Crossspectrum`, this type does not require evenly
sampled data and uses the LSFT instead of the FFT.

Mirrors Python Stingray's `LombScargleCrossspectrum`.

# Fields
- `freq::Vector{T}`: Mid-bin frequencies of the Fourier transform.
- `power::Vector{<:Union{T, Complex{T}}}`: Normalized cross-spectral powers.
- `power_err::Vector{<:Union{T, Complex{T}}}`: Uncertainties on `power` (≈ `power / √m`).
- `unnorm_power::Vector{Complex{T}}`: Unnormalized cross-spectral powers.
- `unnorm_power_err::Vector{Complex{T}}`: Uncertainties on `unnorm_power`.
- `df::T`: Frequency resolution (Hz).
- `dt::T`: Time resolution of the input data (s).
- `n::Int`: Number of data points in the light curves.
- `m::Int`: Number of averaged spectra (always 1 for single-segment LS).
- `k::Int`: Rebinning factor (1 if not rebinned).
- `nphots1::T`: Total number of photons in light curve 1.
- `nphots2::T`: Total number of photons in light curve 2.
- `norm::String`: Normalization: "frac", "abs", "leahy", or "none".
- `power_type::String`: Power representation: "all", "real", or "absolute".
- `fullspec::Bool`: Whether negative frequencies are included.
- `method::String`: LSFT method used: "fast" or "slow".
- `oversampling::Int`: Oversampling factor (for the fast algorithm).
- `variance1::T`: Variance of light curve 1 (count rate if Poisson).
- `variance2::T`: Variance of light curve 2 (count rate if Poisson).
- `err_dist::String`: Error distribution: "poisson" or "gauss".
- `type::String`: Object type identifier, always "crossspectrum".
"""
struct LombScargleCrossspectrum{T<:Real} <: AbstractCrossspectrum
    freq::Vector{T}
    power::Union{Vector{T}, Vector{Complex{T}}}
    power_err::Union{Vector{T}, Vector{Complex{T}}}
    unnorm_power::Vector{Complex{T}}
    unnorm_power_err::Vector{Complex{T}}
    df::T
    dt::T
    n::Int
    m::Int
    k::Int
    nphots1::T
    nphots2::T
    norm::String
    power_type::String
    fullspec::Bool
    method::String
    oversampling::Int
    variance1::T
    variance2::T
    err_dist::String
    type::String
end

function Base.show(io::IO, lscs::LombScargleCrossspectrum{T}) where T
    print(io, "LombScargleCrossspectrum($(lscs.norm), $(lscs.method), ",
          "$(length(lscs.freq)) freq bins, ",
          "df=$(round(lscs.df, sigdigits=4)) Hz, ",
          "freq range=[$(round(lscs.freq[1], sigdigits=4)), ",
          "$(round(lscs.freq[end], sigdigits=4))] Hz)")
end

# ──────────────────────────────────────────────────────────────────────────────
# Internal: LSFT cross-spectrum computation
# ──────────────────────────────────────────────────────────────────────────────

"""
    _ls_cross(lc1, lc2, freq; fullspec=false, method="fast", oversampling=5)

Compute the unnormalized Lomb-Scargle cross spectrum of two light curves.

Returns `(freq, cross)` where `cross = F₁ × conj(F₂)`.

# Arguments
- `lc1::LightCurve`: Light curve from channel 1.
- `lc2::LightCurve`: Light curve from channel 2.
- `freq::AbstractVector{<:Real}`: Frequency grid.

# Keyword Arguments
- `fullspec::Bool=false`: Include negative frequencies via Hermitian symmetry.
- `method::String="fast"`: LSFT method ("fast" or "slow").
- `oversampling::Int=5`: Oversampling factor (fast method only).
"""
function _ls_cross(lc1::LightCurve, lc2::LightCurve,
                   freq::AbstractVector{<:Real};
                   fullspec::Bool=false,
                   method::String="fast",
                   oversampling::Int=5)
    y1 = Float64.(lc1.counts)
    y2 = Float64.(lc2.counts)
    t1 = Float64.(lc1.time)
    t2 = Float64.(lc2.time)

    if method == "slow"
        lsft1 = lsft_slow(y1, t1, freq)
        lsft2 = lsft_slow(y2, t2, freq)
    elseif method == "fast"
        lsft1 = lsft_fast(y1, t1, freq; oversampling=oversampling)
        lsft2 = lsft_fast(y2, t2, freq; oversampling=oversampling)
    else
        throw(ValueError("method must be one of ['fast','slow']"))
    end

    if fullspec
        lsft1, _ = impose_symmetry_lsft(lsft1, sum(y1), Base.length(y1), freq)
        lsft2, freq = impose_symmetry_lsft(lsft2, sum(y2), Base.length(y2), freq)
    end

    cross = lsft1 .* conj.(lsft2)
    return Float64.(freq), cross
end

"""
    _normalize_crossspectrum(unnorm_power, lc1, lc2, norm, n)

Normalize the cross spectrum using `normalize_periodograms` from fourier.jl.

Computes the mean flux as the geometric mean of the two light curves'
mean count rates, matching Python Stingray's normalization approach.
"""
function _normalize_crossspectrum(unnorm_power::AbstractVector{<:Complex},
                                   lc1::LightCurve, lc2::LightCurve,
                                   norm::String, n::Int)
    dt_val = lc1.dt isa AbstractVector ? lc1.dt[1] : lc1.dt
    nphots1 = Float64(sum(lc1.counts))
    nphots2 = Float64(sum(lc2.counts))
    nphots = sqrt(nphots1 * nphots2)

    mean1 = nphots1 / n
    mean2 = nphots2 / n
    mean_flux = sqrt(mean1 * mean2)

    # Determine variance for Leahy normalization
    variance = nothing
    if lc1.err_method == :gaussian
        variance = sqrt(Statistics.var(Float64.(lc1.counts)) *
                        Statistics.var(Float64.(lc2.counts)))
    end

    normalized = normalize_periodograms(
        unnorm_power, Float64(dt_val), n;
        mean_flux=mean_flux, n_ph=nphots,
        norm=norm, variance=variance, power_type="all"
    )

    return normalized
end

"""
    lscrossspectrum_from_lightcurve(lc1, lc2; kwargs...)

Core function that orchestrates the full Lomb-Scargle cross-spectrum pipeline.

1. Builds the frequency grid via `autofrequency`
2. Computes the unnormalized cross spectrum via `_ls_cross`
3. Normalizes the power
4. Applies `power_type` filtering
5. Constructs and returns a `LombScargleCrossspectrum`

Mirrors Python Stingray's `lscrossspectrum_from_lightcurve`.
"""
function lscrossspectrum_from_lightcurve(lc1::LightCurve, lc2::LightCurve;
                                         norm::String="none",
                                         power_type::String="all",
                                         fullspec::Bool=false,
                                         min_freq::Union{Nothing, Real}=nothing,
                                         max_freq::Union{Nothing, Real}=nothing,
                                         df::Union{Nothing, Real}=nothing,
                                         method::String="fast",
                                         oversampling::Int=5,
                                         nyquist_factor::Real=1)
    # Compute time span and dt
    t_length = max(lc1.time[end], lc2.time[end]) - min(lc1.time[1], lc2.time[1])
    dt_val1 = lc1.dt isa AbstractVector ? lc1.dt[1] : lc1.dt
    dt_val2 = lc2.dt isa AbstractVector ? lc2.dt[1] : lc2.dt
    dt_val = min(dt_val1, dt_val2)

    # Build frequency grid
    freq = autofrequency(;
        min_freq=min_freq, max_freq=max_freq, df=df,
        dt=dt_val, length=t_length, nyquist_factor=nyquist_factor
    )

    # Compute unnormalized cross spectrum
    freq, cross = _ls_cross(lc1, lc2, freq;
                             fullspec=fullspec, method=method,
                             oversampling=oversampling)

    n = Base.length(lc1.time)

    # Normalize
    power = _normalize_crossspectrum(cross, lc1, lc2, norm, n)

    # Apply power_type
    if power_type == "real"
        power = real.(power)
        unnorm_display = real.(cross)
    elseif power_type == "absolute"
        power = abs.(power)
        unnorm_display = abs.(cross)
    else
        unnorm_display = cross
    end

    # Compute metadata
    nphots1 = Float64(sum(lc1.counts))
    nphots2 = Float64(sum(lc2.counts))
    freq_df = Base.length(freq) > 1 ? Float64(freq[2] - freq[1]) : 1.0 / t_length

    # Variance
    if lc1.err_method == :gaussian
        variance1 = Float64(Statistics.mean(lc1.count_error))^2
        variance2 = Float64(Statistics.mean(lc2.count_error))^2
        err_dist = "gauss"
    else
        variance1 = Float64(Statistics.mean(Float64.(lc1.counts)))
        variance2 = Float64(Statistics.mean(Float64.(lc2.counts)))
        err_dist = "poisson"
    end

    # Error estimation: power / sqrt(m), m=1 for single-segment LS
    m = 1
    power_err = power ./ sqrt(m)
    unnorm_power_err = cross ./ sqrt(m)

    return LombScargleCrossspectrum{Float64}(
        Float64.(freq),
        power_type == "all" ? Complex{Float64}.(power) :
            (power_type == "real" ? Float64.(power) : Float64.(power)),
        power_type == "all" ? Complex{Float64}.(power_err) :
            (power_type == "real" ? Float64.(power_err) : Float64.(power_err)),
        Complex{Float64}.(cross),
        Complex{Float64}.(unnorm_power_err),
        freq_df, Float64(dt_val), n, m, 1,
        nphots1, nphots2,
        norm, power_type, fullspec,
        method, oversampling,
        variance1, variance2,
        err_dist, "crossspectrum"
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# Public constructors
# ──────────────────────────────────────────────────────────────────────────────

"""
    _validate_ls_inputs(norm, power_type, method, oversampling, fullspec,
                        min_freq, max_freq)

Validate inputs for Lomb-Scargle cross-spectrum / power-spectrum constructors.
Mirrors Python Stingray's `LombScargleCrossspectrum.initial_checks`.
"""
function _validate_ls_inputs(norm::String, power_type::String, method::String,
                              oversampling::Int, fullspec::Bool,
                              min_freq, max_freq)
    if !(norm in ["frac", "abs", "leahy", "none"])
        throw(ValueError("norm must be one of ['frac','abs','leahy','none']"))
    end
    if !(power_type in ["all", "absolute", "real"])
        throw(ValueError("power_type must be one of ['all','absolute','real']"))
    end
    if !(method in ["fast", "slow"])
        throw(ValueError("method must be one of ['fast','slow']"))
    end
    if oversampling < 1
        throw(ArgumentError("oversampling must be ≥ 1"))
    end
    if !isnothing(min_freq) && min_freq < 0
        throw(ValueError("min_freq must be non-negative"))
    end
    if !isnothing(max_freq) && max_freq < 0
        throw(ValueError("max_freq must be non-negative"))
    end
    if !isnothing(max_freq) && !isnothing(min_freq) && max_freq <= min_freq
        throw(ValueError("max_freq must be greater than min_freq"))
    end
end

"""
    LombScargleCrossspectrum(lc1::LightCurve, lc2::LightCurve; kwargs...)

Construct a `LombScargleCrossspectrum` from two `LightCurve` objects.

# Keyword Arguments
- `norm::String="none"`: Normalization ("frac", "abs", "leahy", "none").
- `power_type::String="all"`: Power type ("all", "real", "absolute").
- `fullspec::Bool=false`: Include negative frequencies.
- `min_freq::Union{Nothing, Real}=nothing`: Minimum frequency.
- `max_freq::Union{Nothing, Real}=nothing`: Maximum frequency.
- `df::Union{Nothing, Real}=nothing`: Frequency resolution.
- `method::String="fast"`: LSFT method ("fast" or "slow").
- `oversampling::Int=5`: Oversampling factor.
"""
function LombScargleCrossspectrum(lc1::LightCurve, lc2::LightCurve;
                                   norm::String="none",
                                   power_type::String="all",
                                   fullspec::Bool=false,
                                   min_freq::Union{Nothing, Real}=nothing,
                                   max_freq::Union{Nothing, Real}=nothing,
                                   df::Union{Nothing, Real}=nothing,
                                   method::String="fast",
                                   oversampling::Int=5)
    _validate_ls_inputs(norm, power_type, method, oversampling, fullspec,
                         min_freq, max_freq)

    return lscrossspectrum_from_lightcurve(lc1, lc2;
        norm=norm, power_type=power_type, fullspec=fullspec,
        min_freq=min_freq, max_freq=max_freq, df=df,
        method=method, oversampling=oversampling)
end

"""
    LombScargleCrossspectrum(ev1::EventList, ev2::EventList; dt, kwargs...)

Construct a `LombScargleCrossspectrum` from two `EventList` objects.

Converts event lists to light curves using `create_lightcurve`, then
delegates to the `LightCurve` constructor.

# Arguments
- `ev1::EventList`: Event list for channel 1.
- `ev2::EventList`: Event list for channel 2.

# Keyword Arguments
- `dt::Real`: Time resolution for binning events into light curves. **Required.**
- All other kwargs are passed to the `LightCurve` constructor.
"""
function LombScargleCrossspectrum(ev1::EventList, ev2::EventList;
                                   dt::Real,
                                   norm::String="none",
                                   power_type::String="all",
                                   fullspec::Bool=false,
                                   min_freq::Union{Nothing, Real}=nothing,
                                   max_freq::Union{Nothing, Real}=nothing,
                                   df::Union{Nothing, Real}=nothing,
                                   method::String="fast",
                                   oversampling::Int=5)
    lc1 = create_lightcurve(ev1, dt)
    lc2 = create_lightcurve(ev2, dt)

    return LombScargleCrossspectrum(lc1, lc2;
        norm=norm, power_type=power_type, fullspec=fullspec,
        min_freq=min_freq, max_freq=max_freq, df=df,
        method=method, oversampling=oversampling)
end

# ──────────────────────────────────────────────────────────────────────────────
# LombScarglePowerspectrum
# ──────────────────────────────────────────────────────────────────────────────

"""
    LombScarglePowerspectrum{T<:Real} <: AbstractPowerspectrum

A power spectrum (periodogram) computed from an unevenly sampled time series
using the Lomb-Scargle Fourier Transform.

Mathematically equivalent to a `LombScargleCrossspectrum` of a signal with
itself. The `nphots` field stores the total number of photons (same as
`nphots1` in the cross-spectrum).

Mirrors Python Stingray's `LombScarglePowerspectrum`.

# Fields
Same as `LombScargleCrossspectrum`, plus:
- `nphots::T`: Total number of photons (equals `nphots1`).
"""
struct LombScarglePowerspectrum{T<:Real} <: AbstractPowerspectrum
    freq::Vector{T}
    power::Union{Vector{T}, Vector{Complex{T}}}
    power_err::Union{Vector{T}, Vector{Complex{T}}}
    unnorm_power::Vector{Complex{T}}
    unnorm_power_err::Vector{Complex{T}}
    df::T
    dt::T
    n::Int
    m::Int
    k::Int
    nphots::T
    nphots1::T
    nphots2::T
    norm::String
    power_type::String
    fullspec::Bool
    method::String
    oversampling::Int
    variance1::T
    variance2::T
    err_dist::String
    type::String
end

function Base.show(io::IO, lsps::LombScarglePowerspectrum{T}) where T
    print(io, "LombScarglePowerspectrum($(lsps.norm), $(lsps.method), ",
          "$(length(lsps.freq)) freq bins, ",
          "df=$(round(lsps.df, sigdigits=4)) Hz, ",
          "freq range=[$(round(lsps.freq[1], sigdigits=4)), ",
          "$(round(lsps.freq[end], sigdigits=4))] Hz)")
end

"""
    LombScarglePowerspectrum(lc::LightCurve; kwargs...)

Construct a `LombScarglePowerspectrum` from a `LightCurve` object.

Internally computes the cross-spectrum of the light curve with itself,
then wraps the result as a power spectrum.

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `power_type::String="all"`: Power type ("all", "real", "absolute").
- `fullspec::Bool=false`: Include negative frequencies.
- `min_freq::Union{Nothing, Real}=nothing`: Minimum frequency.
- `max_freq::Union{Nothing, Real}=nothing`: Maximum frequency.
- `df::Union{Nothing, Real}=nothing`: Frequency resolution.
- `method::String="fast"`: LSFT method ("fast" or "slow").
- `oversampling::Int=5`: Oversampling factor.
"""
function LombScarglePowerspectrum(lc::LightCurve;
                                   norm::String="frac",
                                   power_type::String="all",
                                   fullspec::Bool=false,
                                   min_freq::Union{Nothing, Real}=nothing,
                                   max_freq::Union{Nothing, Real}=nothing,
                                   df::Union{Nothing, Real}=nothing,
                                   method::String="fast",
                                   oversampling::Int=5)
    _validate_ls_inputs(norm, power_type, method, oversampling, fullspec,
                         min_freq, max_freq)

    # Compute as cross-spectrum of signal with itself
    cs = lscrossspectrum_from_lightcurve(lc, lc;
        norm=norm, power_type=power_type, fullspec=fullspec,
        min_freq=min_freq, max_freq=max_freq, df=df,
        method=method, oversampling=oversampling)

    return LombScarglePowerspectrum{Float64}(
        cs.freq, cs.power, cs.power_err,
        cs.unnorm_power, cs.unnorm_power_err,
        cs.df, cs.dt, cs.n, cs.m, cs.k,
        cs.nphots1,  # nphots = nphots1 for auto-spectrum
        cs.nphots1, cs.nphots2,
        cs.norm, cs.power_type, cs.fullspec,
        cs.method, cs.oversampling,
        cs.variance1, cs.variance2,
        cs.err_dist, "powerspectrum"
    )
end

"""
    LombScarglePowerspectrum(ev::EventList; dt, kwargs...)

Construct a `LombScarglePowerspectrum` from an `EventList` object.

Converts the event list to a light curve, then delegates to the
`LightCurve` constructor.

# Keyword Arguments
- `dt::Real`: Time resolution for binning. **Required.**
- All other kwargs are passed to the `LightCurve` constructor.
"""
function LombScarglePowerspectrum(ev::EventList;
                                   dt::Real,
                                   norm::String="frac",
                                   power_type::String="all",
                                   fullspec::Bool=false,
                                   min_freq::Union{Nothing, Real}=nothing,
                                   max_freq::Union{Nothing, Real}=nothing,
                                   df::Union{Nothing, Real}=nothing,
                                   method::String="fast",
                                   oversampling::Int=5)
    lc = create_lightcurve(ev, dt)

    return LombScarglePowerspectrum(lc;
        norm=norm, power_type=power_type, fullspec=fullspec,
        min_freq=min_freq, max_freq=max_freq, df=df,
        method=method, oversampling=oversampling)
end
