module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames, FFTW, NaNMath, FITSIO, Intervals, LinearAlgebra, Colors, Plots, StatsPlots
using ProgressBars: tqdm as show_progress

include("fourier.jl")
include("gti.jl")
include("utils.jl")
include("power_colors.jl")

# --- Export from fourier.jl ---
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

# --- Export from gti.jl ---
export load_gtis
export get_total_gti_length
export create_gti_mask
export create_gti_from_condition
export operations_on_gtis
export get_btis
export time_intervals_from_gtis
export bin_intervals_from_gtis

# --- Export from power_colors.jl ---
export plot_power_colors
export hue_from_power_color
export plot_hues
export integrate_power_in_frequency_range
export power_color
export hue_from_logpower_color


include("events.jl")
export readevents, EventList, DictMetadata

end
