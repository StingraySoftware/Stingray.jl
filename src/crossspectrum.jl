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
- `m::Int`: Number of averaged cross-spectra.
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
    print(io, "Crossspectrum($(cs.norm), $(cs.m) segments, ",
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

    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    # Extract spectral data from the DataFrame
    freq = Float64.(result[!, "freq"])
    power = Complex{Float64}.(result[!, "power"])
    unnorm_power = Complex{Float64}.(result[!, "unnorm_power"])

    # Compute metadata from input parameters
    n_bin = round(Int, segment_size / dt)
    dt_adj = segment_size / n_bin
    df_val = 1.0 / segment_size
    m_val = _count_valid_segments(common_gti, segment_size)

    power_err = power ./ sqrt(m_val)
    unnorm_power_err = unnorm_power ./ sqrt(m_val)

    # Extract auxiliary PDS if computed
    _pds1 = "pds1" in names(result) ? Float64.(real.(result[!, "pds1"])) : nothing
    _pds2 = "pds2" in names(result) ? Float64.(real.(result[!, "pds2"])) : nothing

    # Compute photon statistics
    nphots1_val = Float64(length(times(ev1))) / m_val
    nphots2_val = Float64(length(times(ev2))) / m_val
    nphots_val = sqrt(nphots1_val * nphots2_val)
    countrate1_val = nphots1_val / segment_size
    countrate2_val = nphots2_val / segment_size

    return Crossspectrum{Float64}(
        freq, power, power_err, unnorm_power, unnorm_power_err,
        _pds1, _pds2,
        df_val, dt_adj, n_bin, m_val,
        1,  # k = 1 (not rebinned)
        nphots_val, nphots1_val, nphots2_val,
        countrate1_val, countrate2_val,
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

    if isnothing(result)
        throw(ArgumentError("No valid segments found. Check GTIs and segment_size."))
    end

    freq = Float64.(result[!, "freq"])
    power = Complex{Float64}.(result[!, "power"])
    unnorm_power = Complex{Float64}.(result[!, "unnorm_power"])

    n_bin = round(Int, segment_size / dt_val)
    dt_adj = segment_size / n_bin
    df_val = 1.0 / segment_size
    m_val = _count_valid_segments(common_gti, segment_size)

    power_err = power ./ sqrt(m_val)
    unnorm_power_err = unnorm_power ./ sqrt(m_val)

    _pds1 = "pds1" in names(result) ? Float64.(real.(result[!, "pds1"])) : nothing
    _pds2 = "pds2" in names(result) ? Float64.(real.(result[!, "pds2"])) : nothing

    nphots1_val = Float64(sum(lc1.counts)) / m_val
    nphots2_val = Float64(sum(lc2.counts)) / m_val
    nphots_val = sqrt(nphots1_val * nphots2_val)
    countrate1_val = nphots1_val / segment_size
    countrate2_val = nphots2_val / segment_size

    # Determine error distribution from LightCurve err_method
    err_dist = lc1.err_method == :gaussian ? "gauss" : "poisson"

    return Crossspectrum{Float64}(
        freq, power, power_err, unnorm_power, unnorm_power_err,
        _pds1, _pds2,
        df_val, dt_adj, n_bin, m_val,
        1,
        nphots_val, nphots1_val, nphots2_val,
        countrate1_val, countrate2_val,
        norm, power_type, fullspec,
        common_gti, Float64(segment_size),
        err_dist, "crossspectrum"
    )
end
