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
