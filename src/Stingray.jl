module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames, FFTW, Metadata, NaNMath
using ProgressBars: tqdm as show_progress 

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
export time_intervals_from_gtis
export bin_intervals_from_gtis

include("utils.jl")

end 
