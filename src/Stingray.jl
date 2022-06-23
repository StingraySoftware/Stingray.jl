module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames, FFTW, Metadata
using ProgressBars: tqdm as show_progress 

include("fourier.jl")
export positive_fft_bins
export normalize_periodograms
export normalize_abs
export normalize_frac
export poisson_level
export bias_term
export raw_coherence
export estimate_intrinsic_coherence
export get_average_ctrate
export get_flux_iterable_from_segments
export avg_pds_from_iterable
export avg_pds_from_events
export avg_cs_from_events
export normalize_leahy_from_variance

include("gti.jl")
export generate_indices_of_segment_boundaries_binned
export generate_indices_of_segment_boundaries_unbinned
export bin_intervals_from_gtis
export calculate_segment_bin_start
export time_intervals_from_gtis

include("utils.jl")
export sum_if_not_none_or_initialize

end 
