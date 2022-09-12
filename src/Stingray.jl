module Stingray

using ResumableFunctions, StatsBase, Statistics, DataFrames
using FFTW, Metadata, NaNMath, FITSIO, Intervals, Parameters
using ProgressBars: tqdm as show_progress 

include("fourier.jl")
export positive_fft_bins
export poisson_level
export normalize_periodograms
export raw_coherence
export estimate_intrinsic_coherence
export error_on_averaged_cross_spectrum
export get_average_ctrate
export avg_pds_from_events
export avg_cs_from_events

include("gti.jl")
export load_gtis
export get_total_gti_length
export create_gti_mask
export create_gti_from_condition
export operations_on_gtis
export get_btis

include("lightcurve.jl")
export LightCurve
export make_lightcurves
export rebin

include("events.jl")
export EventList
export to_lc
export from_lc

include("io.jl")
export load_events_from_fits

include("utils.jl")

end 
