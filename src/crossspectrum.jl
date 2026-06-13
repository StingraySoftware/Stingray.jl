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
