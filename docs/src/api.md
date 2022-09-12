# API Reference

## EventList Operations

```@docs
EventList
read
write
from_lc
to_lc
join
sort
sort!
Stingray.apply_mask
```
## LightCurve Operations
```@docs
LightCurve
make_lightcurves
rebin
```

## GTI Operations
```@docs
load_gtis
get_total_gti_length
create_gti_mask
create_gti_from_condition
operations_on_gtis
get_btis
```

## Fourier Operations
```@docs
positive_fft_bins
poisson_level
normalize_periodograms
raw_coherence
estimate_intrinsic_coherence
error_on_averaged_cross_spectrum
get_average_ctrate
avg_pds_from_events
avg_cs_from_events
```
