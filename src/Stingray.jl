module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames
using FFTW, NaNMath, FITSIO, Intervals
using ProgressBars: tqdm as show_progress
using DocStringExtensions
using LinearAlgebra
using Random
using RecipesBase

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
    filter_on!,
    read_gti_from_fits,
    gti_info,
    gti_exposure,
    gti,
    has_gti,
    extract_timing_keywords

include("missionSupport.jl")
export MissionSupport,
    get_mission_support,
    apply_calibration,
    patch_mission_info,
    SIMPLE_CALIBRATION_FUNCS,
    interpret_fits_data!,
    AbstractMissionSupport

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
include("utils.jl")

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
export avg_pds_from_events, avg_pds_from_iterable
export avg_cs_from_events, avg_cs_from_iterables, avg_cs_from_iterables_quick
export avg_pds_from_eventlist,
    avg_cs_from_eventlists, avg_pds_from_lightcurve, avg_cs_from_lightcurves
export get_norm_label, get_poisson_level, extract_gti

include("gti.jl")
export load_gtis
export get_total_gti_length
export create_gti_mask
export create_gti_from_condition
export operations_on_gtis
export get_btis
export time_intervals_from_gtis
export bin_intervals_from_gtis
export apply_gtis
export fill_bad_time_intervals!
export create_filtered_lightcurve
export check_gtis
export split_by_gtis, intersect_gtis,get_gti_lengths
include("powerspectrum.jl")
export Powerspectrum
export AveragedPowerspectrum, AbstractPowerSpectrum
include("crossspectrum.jl")
export CrossSpectrum, AveragedCrossSpectrum
export is_averaged, is_single, theoretical_noise_level, fill_errors!
export white_noise_level, noise_corrected_power, signal_to_noise_ratio
export detect_aliasing, coherence, phase_lag, time_lag, noise_properties
export significant_frequencies, get_noise_level, quality_metrics
export rebin, rebin_log, geometric_rebin, adaptive_rebin
export is_rebinned, effective_samples_per_bin, AbstractCrossSpectrum
include("plotting/plots_recipes_lightcurve.jl")
export create_segments
include("plotting/plots_recipes_gti.jl")
export BTIAnalysisPlot
include("plotting/plots_recipes_crossspectrum.jl")
include("plotting/plots_recipes_powerspectrum.jl")

end
