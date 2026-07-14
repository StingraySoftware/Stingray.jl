"""
Abstract base type for all cross-spectral objects.

All cross-spectral types (Crossspectrum, AveragedCrossspectrum, etc.)
should subtype this. Mirrors the Python Stingray class hierarchy where
`Crossspectrum` is the base class.
"""
abstract type AbstractCrossspectrum end

"""
Abstract base type for all power spectral objects.

Subtypes `AbstractCrossspectrum` because a power spectrum is mathematically
a cross-spectrum of a signal with itself. This mirrors the Python Stingray
hierarchy where `Powerspectrum(Crossspectrum)`.
"""
abstract type AbstractPowerspectrum <: AbstractCrossspectrum end

"""
    Crossspectrum{T<:Real} <: AbstractCrossspectrum

A cross spectrum computed from two time series (event lists or light curves).

The cross spectrum is the Fourier-domain equivalent of the cross-correlation
function and encodes information about the coherence and time lags between
two signals at each Fourier frequency.

# Fields
- `freq::Vector{T}`: Mid-bin frequencies of the Fourier transform.
- `power::Vector{Complex{T}}`: Normalized cross-spectral powers (complex).
- `power_err::Vector{Complex{T}}`: Uncertainties on `power`, approximated as `power / √m`.
- `unnorm_power::Vector{Complex{T}}`: Unnormalized cross-spectral powers.
- `unnorm_power_err::Vector{Complex{T}}`: Uncertainties on `unnorm_power`.
- `pds1::Union{Nothing, Vector{T}}`: Auxiliary power spectrum of signal 1 (if computed).
- `pds2::Union{Nothing, Vector{T}}`: Auxiliary power spectrum of signal 2 (if computed).
- `df::T`: Frequency resolution (Hz).
- `dt::T`: Time resolution of the input data (s).
- `n::Int`: Number of bins per segment.
- `m::Union{Int, Vector{Int}}`: Number of averaged cross-spectra.
- `k::Union{Int, Vector{Int}}`: Rebinning factor (1 if not rebinned, vector after log rebinning).
- `nphots::T`: Geometric mean of total photons: √(nphots1 * nphots2).
- `nphots1::T`: Total number of photons in signal 1.
- `nphots2::T`: Total number of photons in signal 2.
- `countrate1::T`: Mean count rate of signal 1 (counts/s).
- `countrate2::T`: Mean count rate of signal 2 (counts/s).
- `norm::String`: Normalization applied: "frac", "abs", "leahy", or "none".
- `power_type::String`: Type of power stored: "all" (complex), "real", or "abs".
- `fullspec::Bool`: Whether negative frequencies are included.
- `gti::Matrix{T}`: Good Time Intervals, Nx2 matrix of [start, stop] pairs.
- `segment_size::T`: Length of each segment in seconds.
- `err_dist::String`: Error distribution assumption: "poisson" or "gauss".
- `type::String`: Object type identifier, always "crossspectrum".
"""
struct Crossspectrum{T<:Real} <: AbstractCrossspectrum
    freq::Vector{T}
    power::Vector{Complex{T}}
    power_err::Vector{Complex{T}}
    unnorm_power::Vector{Complex{T}}
    unnorm_power_err::Vector{Complex{T}}
    pds1::Union{Nothing, Vector{T}}
    pds2::Union{Nothing, Vector{T}}
    df::T
    dt::T
    n::Int
    m::Int
    k::Union{Int, Vector{Int}}
    nphots::T
    nphots1::T
    nphots2::T
    countrate1::T
    countrate2::T
    norm::String
    power_type::String
    fullspec::Bool
    gti::Matrix{T}
    segment_size::T
    err_dist::String
    type::String
end

function Base.show(io::IO, cs::Crossspectrum{T}) where T
    m_display = cs.m isa Int ? cs.m : "$(length(cs.m))-element"
    print(io, "Crossspectrum($(cs.norm), $(m_display) segments, ",
          "$(length(cs.freq)) freq bins, ",
          "df=$(round(cs.df, sigdigits=4)) Hz, ",
          "freq range=[$(round(cs.freq[1], sigdigits=4)), ",
          "$(round(cs.freq[end], sigdigits=4))] Hz)")
end

"""
    _count_valid_segments(gti::Matrix, segment_size::Real) -> Int

Count the number of valid segments that fit within the GTI intervals.
Mirrors the segment counting logic in `get_flux_iterable_from_segments`.
"""
function _count_valid_segments(gti::Matrix, segment_size::Real)
    m = 0
    for i in 1:size(gti, 1)
        m += floor(Int, (gti[i, 2] - gti[i, 1]) / segment_size)
    end
    return max(m, 1)
end

"""
    _compute_spectral_params(segment_size, dt, gti) -> NamedTuple

Compute shared spectral metadata (n_bin, adjusted dt, df, m) from segment
parameters and GTIs. Used internally by `Crossspectrum` and `Powerspectrum`
constructors.
"""
function _compute_spectral_params(segment_size::Real, dt::Real, gti::Matrix)
    n_bin = round(Int, segment_size / dt)
    dt_adj = segment_size / n_bin
    df = 1.0 / segment_size
    m = _count_valid_segments(gti, segment_size)
    return (; n_bin, dt=dt_adj, df, m)
end

"""
    _extract_spectral_result(result; complex_power=true) -> NamedTuple

Extract freq, power, unnorm_power, and auxiliary PDS columns from a DataFrame
returned by `avg_cs_from_events` or `avg_pds_from_events`. Throws if `result`
is `nothing`.  Set `complex_power=false` for real-valued power (Powerspectrum).
"""
function _extract_spectral_result(result; complex_power::Bool=true)
    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    freq = Float64.(result[!, "freq"])
    if complex_power
        power = Complex{Float64}.(result[!, "power"])
        unnorm_power = Complex{Float64}.(result[!, "unnorm_power"])
    else
        power = Float64.(real.(result[!, "power"]))
        unnorm_power = Float64.(real.(result[!, "unnorm_power"]))
    end

    pds1 = "pds1" in names(result) ? Float64.(real.(result[!, "pds1"])) : nothing
    pds2 = "pds2" in names(result) ? Float64.(real.(result[!, "pds2"])) : nothing

    return (; freq, power, unnorm_power, pds1, pds2)
end

"""
    _compute_power_errors(power, unnorm_power, m) -> (power_err, unnorm_power_err)

Estimate spectral uncertainties as `power / √m`. Supports both scalar `m`
(uniform averaging) and vector `m` (variable averaging after log rebinning).
"""
_compute_power_errors(power, unnorm_power, m::Int) =
    (power ./ sqrt(m), unnorm_power ./ sqrt(m))
_compute_power_errors(power, unnorm_power, m::Vector{Int}) =
    (power ./ sqrt.(m), unnorm_power ./ sqrt.(m))

"""
    _compute_photon_stats(nphots1, nphots2, m, segment_size) -> NamedTuple

Compute per-segment photon statistics for a cross-spectrum (two signals).
"""
function _compute_photon_stats(nphots1::Real, nphots2::Real, m::Int, segment_size::Real)
    n1 = Float64(nphots1) / m
    n2 = Float64(nphots2) / m
    return (; nphots=sqrt(n1 * n2), nphots1=n1, nphots2=n2,
              countrate1=n1 / segment_size, countrate2=n2 / segment_size)
end

"""
    _compute_photon_stats(nphots_total, m, segment_size) -> NamedTuple

Compute per-segment photon statistics for a power spectrum (single signal).
"""
function _compute_photon_stats(nphots_total::Real, m::Int, segment_size::Real)
    return (; nphots=Float64(nphots_total) / m)
end

"""
    Crossspectrum(ev1::EventList, ev2::EventList, segment_size::Real, dt::Real; kwargs...)

Construct a `Crossspectrum` from two `EventList` objects.

Calls `avg_cs_from_events` and computes metadata from input parameters.
Error bars are computed as `power / √m`.

# Arguments
- `ev1::EventList`: Event list for the subject band.
- `ev2::EventList`: Event list for the reference band.
- `segment_size::Real`: Length of each segment in seconds.
- `dt::Real`: Time resolution (sets Nyquist frequency at 0.5/dt Hz).

# Keyword Arguments
- `norm::String="frac"`: Normalization ("frac", "abs", "leahy", "none").
- `use_common_mean::Bool=true`: Use the mean from the full light curve.
- `silent::Bool=false`: Suppress progress bars.
- `fullspec::Bool=false`: Include negative frequencies.
- `power_type::String="all"`: Power type ("all", "real", "abs").
- `return_auxil::Bool=false`: Compute auxiliary PDS for each signal.
"""
function Crossspectrum(ev1::EventList, ev2::EventList,
                       segment_size::Real, dt::Real;
                       norm::String="frac",
                       use_common_mean::Bool=true,
                       silent::Bool=false,
                       fullspec::Bool=false,
                       power_type::String="all",
                       return_auxil::Bool=false)

    common_gti = has_gti(ev1) && has_gti(ev2) ?
        operations_on_gtis([gti(ev1), gti(ev2)], intersect) :
        (has_gti(ev1) ? gti(ev1) : gti(ev2))

    result = avg_cs_from_events(
        times(ev1), times(ev2), common_gti,
        segment_size, dt;
        norm=norm,
        use_common_mean=use_common_mean,
        silent=silent,
        fullspec=fullspec,
        power_type=power_type,
        return_auxil=return_auxil
    )

    spectral = _extract_spectral_result(result; complex_power=true)
    params = _compute_spectral_params(segment_size, dt, common_gti)
    power_err, unnorm_power_err = _compute_power_errors(
        spectral.power, spectral.unnorm_power, params.m)
    stats = _compute_photon_stats(
        length(times(ev1)), length(times(ev2)), params.m, segment_size)

    return Crossspectrum{Float64}(
        spectral.freq, spectral.power, power_err,
        spectral.unnorm_power, unnorm_power_err,
        spectral.pds1, spectral.pds2,
        params.df, params.dt, params.n_bin, params.m,
        1,  # k = 1 (not rebinned)
        stats.nphots, stats.nphots1, stats.nphots2,
        stats.countrate1, stats.countrate2,
        norm, power_type, fullspec,
        common_gti, Float64(segment_size),
        "poisson", "crossspectrum"
    )
end

"""
    Crossspectrum(lc1::LightCurve, lc2::LightCurve, segment_size::Real; kwargs...)

Construct a `Crossspectrum` from two `LightCurve` objects.

Extracts times, counts, and derives GTIs from the light curves' time ranges.
The dt is taken from `lc.dt` (the bin size stored on the LightCurve struct).

# Arguments
- `lc1::LightCurve`: Light curve for the subject band.
- `lc2::LightCurve`: Light curve for the reference band.
- `segment_size::Real`: Length of each segment in seconds.

# Keyword Arguments
Same as the EventList constructor, except `dt` is derived from the light curve.
"""
function Crossspectrum(lc1::LightCurve, lc2::LightCurve,
                       segment_size::Real;
                       norm::String="frac",
                       use_common_mean::Bool=true,
                       silent::Bool=false,
                       fullspec::Bool=false,
                       power_type::String="all",
                       return_auxil::Bool=false)

    # dt is a direct field on LightCurve, handle Union{T,Vector{T}}
    dt_val = lc1.dt isa AbstractVector ? lc1.dt[1] : lc1.dt

    # Derive GTI from light curve time ranges (LightCurveMetadata has no gti field)
    lc1_gti = [lc1.metadata.time_range[1] lc1.metadata.time_range[2]]
    lc2_gti = [lc2.metadata.time_range[1] lc2.metadata.time_range[2]]
    common_gti = operations_on_gtis([lc1_gti, lc2_gti], intersect)

    result = avg_cs_from_events(
        lc1.time, lc2.time, common_gti,
        segment_size, dt_val;
        norm=norm,
        use_common_mean=use_common_mean,
        silent=silent,
        fullspec=fullspec,
        power_type=power_type,
        fluxes1=Float64.(lc1.counts),
        fluxes2=Float64.(lc2.counts),
        return_auxil=return_auxil
    )

    spectral = _extract_spectral_result(result; complex_power=true)
    params = _compute_spectral_params(segment_size, dt_val, common_gti)
    power_err, unnorm_power_err = _compute_power_errors(
        spectral.power, spectral.unnorm_power, params.m)
    stats = _compute_photon_stats(
        sum(lc1.counts), sum(lc2.counts), params.m, segment_size)

    err_dist = lc1.err_method == :gaussian ? "gauss" : "poisson"

    return Crossspectrum{Float64}(
        spectral.freq, spectral.power, power_err,
        spectral.unnorm_power, unnorm_power_err,
        spectral.pds1, spectral.pds2,
        params.df, params.dt, params.n_bin, params.m,
        1,
        stats.nphots, stats.nphots1, stats.nphots2,
        stats.countrate1, stats.countrate2,
        norm, power_type, fullspec,
        common_gti, Float64(segment_size),
        err_dist, "crossspectrum"
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# Rebinning utilities
# ──────────────────────────────────────────────────────────────────────────────


"""
    rebin_data(x, y, dx_new; yerr=nothing, method=:mean, dx=nothing)

Rebin data to a new resolution `dx_new`. Either sum or average the data points
in the new bins. Handles partial bins at edges using fractional overlap.

Mirrors Python Stingray's `stingray.utils.rebin_data`.

# Arguments
- `x::AbstractVector`: The dependent variable (e.g., frequencies).
- `y::AbstractVector`: The independent variable to be binned (real or complex).
- `dx_new::Real`: The new resolution of `x`.

# Keyword Arguments
- `yerr::Union{Nothing, AbstractVector}`: Uncertainties on `y`. Propagated via
  quadrature during binning.
- `method::Symbol`: `:mean` (average) or `:sum`. Default `:mean`.
- `dx::Union{Nothing, Real}`: The old resolution. If `nothing`, computed from
  `diff(x)`.

# Returns
- `(xbin, ybin, ybinerr, step_size)` — rebinned x midpoints, rebinned y values,
  propagated errors, and the number of old bins per new bin.
"""
function rebin_data(x::AbstractVector, y::AbstractVector, dx_new::Real;
                    yerr::Union{Nothing, AbstractVector}=nothing,
                    method::Symbol=:mean,
                    dx::Union{Nothing, Real}=nothing)

    y = collect(y)
    if isnothing(yerr)
        yerr_arr = zeros(Float64, length(y))
    else
        yerr_arr = collect(yerr)
    end

    # Determine old resolution
    if isnothing(dx) || dx == 0
        dx_old = diff(x)
    else
        dx_old = [float(dx)]
    end

    if any(dx_new .< dx_old)
        throw(ArgumentError(
            "New frequency resolution must be larger than old frequency resolution."))
    end

    # Left and right bin edges — assumes x values are left bin edges
    xedges = vcat(collect(Float64, x), [Float64(x[end]) + dx_old[end]])

    # New regularly-binned edges
    # Match Python's np.arange(xedges[0], xedges[-1] + dx_new, dx_new)
    xbin_edges = collect(xedges[1]:dx_new:(xedges[end] + dx_new))

    n_new = length(xbin_edges) - 1
    T_out = eltype(y)

    output = zeros(T_out, n_new)
    outputvar = zeros(Float64, n_new)
    step_size = zeros(Float64, n_new)

    # searchsortedfirst matches Python's np.searchsorted (side='left')
    all_x = [searchsortedfirst(xedges, xb) for xb in xbin_edges]

    for i in 1:n_new
        min_ind = all_x[i]
        max_ind = all_x[i + 1]
        xmin = xbin_edges[i]
        xmax = xbin_edges[i + 1]

        # Sum full bins: y[min_ind : max_ind-2] in Julia (matches Python y[min_ind:max_ind-1])
        for j in min_ind:(max_ind - 2)
            if 1 <= j <= length(y)
                output[i] += y[j]
                outputvar[i] += abs(yerr_arr[j])^2
                step_size[i] += 1.0
            end
        end

        # Fractional contribution from the bin straddling the left edge
        if min_ind >= 2 && min_ind - 1 <= length(y)
            prev_dx = xedges[min_ind] - xedges[min_ind - 1]
            prev_frac = (xedges[min_ind] - xmin) / prev_dx
            output[i] += y[min_ind - 1] * prev_frac
            outputvar[i] += (abs(yerr_arr[min_ind - 1]) * prev_frac)^2
            step_size[i] += prev_frac
        end

        # Fractional contribution from the bin straddling the right edge
        if max_ind <= length(xedges) && max_ind - 1 >= 1 && max_ind - 1 <= length(y)
            dx_post = xedges[max_ind] - xedges[max_ind - 1]
            post_frac = (xmax - xedges[max_ind - 1]) / dx_post
            output[i] += y[max_ind - 1] * post_frac
            outputvar[i] += (abs(yerr_arr[max_ind - 1]) * post_frac)^2
            step_size[i] += post_frac
        end
    end

    # Apply method
    if method in (:mean, :avg, :average)
        ybin = output ./ step_size
        ybinerr = sqrt.(outputvar) ./ step_size
    elseif method == :sum
        ybin = output
        ybinerr = sqrt.(outputvar)
    else
        throw(ArgumentError(
            "Method not recognized. Please use :sum or :mean."))
    end

    # Handle non-evenly-divisible segments (matches Python trim logic)
    tseg = length(dx_old) == 1 ?
        Float64(x[end]) - Float64(x[1]) + dx_old[1] :
        Float64(x[end]) - Float64(x[1]) + dx_old[end]

    if (tseg / dx_new) % 1 > 0
        ybin = ybin[1:end-1]
        ybinerr = ybinerr[1:end-1]
        step_size = step_size[1:end-1]
    end

    # Also trim any trailing bins with zero samples (floating-point edge artifacts)
    while length(step_size) > 0 && step_size[end] <= 0
        ybin = ybin[1:end-1]
        ybinerr = ybinerr[1:end-1]
        step_size = step_size[1:end-1]
    end

    # Compute bin midpoints
    n_out = length(ybin)
    xbin = [xbin_edges[i] + dx_new / 2 for i in 1:n_out]

    return xbin, ybin, ybinerr, step_size
end

"""
    _root_squared_mean(x)

Compute √(mean(x²)) — the root-mean-square statistic used for error
propagation during logarithmic rebinning.
"""
_root_squared_mean(x) = sqrt(mean(x .^ 2))

"""
    rebin_data_log(x, y, f; y_err=nothing, dx=nothing)

Logarithmic re-bin of data. Each new bin width grows as `dν_j = dν_{j-1} * (1+f)`.

Mirrors Python Stingray's `stingray.utils.rebin_data_log`.

# Arguments
- `x::AbstractVector`: The dependent variable (e.g., frequencies).
- `y::AbstractVector`: The independent variable to be binned (real or complex).
- `f::Real`: The factor of increase of each bin width relative to the previous one.

# Keyword Arguments
- `y_err::Union{Nothing, AbstractVector}`: Uncertainties on `y`.
- `dx::Union{Nothing, Real}`: Initial bin width. If `nothing`, computed as `median(diff(x))`.

# Returns
- `(xbin, ybin, ybin_err, nsamples)` — rebinned midpoints, rebinned values,
  propagated errors, and number of original samples per bin.
"""
function rebin_data_log(x::AbstractVector, y::AbstractVector, f::Real;
                        y_err::Union{Nothing, AbstractVector}=nothing,
                        dx::Union{Nothing, Real}=nothing)

    x = Float64.(collect(x))
    y = collect(y)
    y_err_arr = isnothing(y_err) ? zeros(Float64, length(y)) : collect(y_err)

    if length(x) != length(y)
        throw(ArgumentError("x and y must be of the same length!"))
    end
    if length(y) != length(y_err_arr)
        throw(ArgumentError("y and y_err must be of the same length!"))
    end

    dx_init = isnothing(dx) ? median(diff(x)) : Float64(dx)

    # Build logarithmically-growing bin edges
    minx = x[1] * 0.5    # frequency to start from
    maxx = x[end]         # maximum frequency to end
    binx = [minx, minx + dx_init]
    dx_curr = dx_init

    while binx[end] <= maxx
        push!(binx, binx[end] + dx_curr * (1.0 + f))
        dx_curr = binx[end] - binx[end-1]
    end

    binx_edges = Float64.(binx)
    n_bins = length(binx_edges) - 1

    # Assign each x value to a bin
    bin_indices = zeros(Int, length(x))
    for (j, xv) in enumerate(x)
        idx = searchsortedlast(binx_edges, xv)
        if idx >= 1 && idx <= n_bins
            bin_indices[j] = idx
        end
    end

    # Compute binned statistics
    xbin = Vector{Float64}()
    ybin_real = Vector{Float64}()
    ybin_imag = Vector{Float64}()
    ybin_err_real = Vector{Float64}()
    ybin_err_imag = Vector{Float64}()
    nsamples = Vector{Int}()

    is_complex = eltype(y) <: Complex

    for i in 1:n_bins
        mask = bin_indices .== i
        count = sum(mask)
        count == 0 && continue

        # Mean of x values in this bin
        push!(xbin, mean(x[mask]))

        # Mean of y values (handle complex)
        vals = y[mask]
        push!(ybin_real, mean(real.(vals)))
        if is_complex
            push!(ybin_imag, mean(imag.(vals)))
        end

        # RMS error propagation
        errs = y_err_arr[mask]
        push!(ybin_err_real, _root_squared_mean(real.(errs)))
        if is_complex
            push!(ybin_err_imag, _root_squared_mean(imag.(errs)))
        end

        push!(nsamples, count)
    end

    # Reconstruct output
    if is_complex
        ybin_out = ybin_real .+ im .* ybin_imag
        ybin_err_out = ybin_err_real .+ im .* ybin_err_imag
    else
        ybin_out = ybin_real
        ybin_err_out = ybin_err_real
    end

    return xbin, ybin_out, ybin_err_out, nsamples
end


# ──────────────────────────────────────────────────────────────────────────────
# Crossspectrum rebinning methods
# ──────────────────────────────────────────────────────────────────────────────

"""
    rebin(cs::Crossspectrum{T}, df_new::Real; f=nothing, method=:mean)

Rebin the cross spectrum to a new frequency resolution `df_new`.

Mirrors Python Stingray's `Crossspectrum.rebin()`.

# Arguments
- `cs::Crossspectrum{T}`: The cross spectrum to rebin.
- `df_new::Real`: The new frequency resolution.

# Keyword Arguments
- `f::Union{Nothing, Real}`: Rebin factor. If specified, `df_new = f * cs.df`.
- `method::Symbol`: `:mean` or `:sum`. Default `:mean`.

# Returns
A new `Crossspectrum` with updated `freq`, `power`, `df`, `m`, and `k`.
"""
function rebin(cs::Crossspectrum{T}, df_new::Real;
               f::Union{Nothing, Real}=nothing,
               method::Symbol=:mean) where T

    if !isnothing(f)
        df_new = f * cs.df
    end
