module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames
using FFTW, NaNMath, FITSIO, Intervals
using ProgressBars: tqdm as show_progress
using DocStringExtensions
using LinearAlgebra

include("fourier.jl")
export positive_fft_bins
export poisson_level
export normalize_abs
export normalize_frac
export normalize_leahy_from_variance
export normalize_periodograms
export bias_term
export raw_coherence
export estimate_intrinsic_coherence
export error_on_averaged_cross_spectrum
export get_average_ctrate
export get_flux_iterable_from_segments
export avg_pds_from_events
export avg_cs_from_events

include("gti.jl")
export load_gtis
export get_total_gti_length
export create_gti_mask
export create_gti_from_condition
export operations_on_gtis
export get_btis
export time_intervals_from_gtis
export bin_intervals_from_gtis

include("utils.jl")

include("events.jl")
export FITSMetadata,
    EventList,
    times,
    energies,
    has_energies,
    filter_time!,
    filter_energy!,
    filter_time,
    filter_energy,
    colnames,
    read_energy_column,
    readevents,
    summary,
    filter_on!

include("lightcurve.jl")
export AbstractLightCurve,
       EventProperty,
       LightCurveMetadata,
       LightCurve,
       calculate_errors,
       set_errors!,
       calculate_errors!,
       create_time_bins,
       bin_events,
       apply_filters,
       calculate_event_properties,
       extract_metadata,
       create_lightcurve,
       rebin

end
