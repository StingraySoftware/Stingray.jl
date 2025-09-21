# Cross-Spectrum

## Overview

Cross-spectrum analysis examines the relationship between two time series in the frequency domain, revealing correlations, phase relationships, and spectral properties that are invisible in individual power spectra. This implementation provides comprehensive tools for creating, analyzing, and visualizing cross-power spectra from astronomical timing data.

### What You'll Learn
- How to create single-segment and averaged cross-spectra
- Understanding different normalization schemes and their applications
- Noise analysis and statistical significance testing
- Rebinning strategies for different scientific goals
- Comprehensive visualization techniques
- Physical interpretation of coherence, phase lags, and time lags

## Core Data Structure

### `CrossSpectrum{T}`

The unified data structure handles both single-segment and averaged cross-spectra:

```julia
mutable struct CrossSpectrum{T} <: AbstractCrossSpectrum{T}
    freq::Vector{T}                    # Frequencies in Hz
    power::Vector{Complex{T}}          # Complex cross power values  
    power_err::Union{Nothing,Vector{T}} # Power errors
    ps1::Vector{T}                     # Auto power spectrum of first signal
    ps2::Vector{T}                     # Auto power spectrum of second signal
    norm::String                       # Normalization type
    power_type::String                 # Power type (all, real, absolute)
    df::T                             # Frequency resolution (Hz)
    nphots1::T                        # Total photons in first signal
    nphots2::T                        # Total photons in second signal
    m::Int                            # Number of segments (1=single, >1=averaged)
    n::Int                            # Number of frequencies
    k::Union{Int,Vector{Int}}         # Rebinning factor per frequency bin
    metadata1::Union{LightCurveMetadata,FITSMetadata}
    metadata2::Union{LightCurveMetadata,FITSMetadata}
    fullspec::Bool                    # Include negative frequencies
    channels_overlap::Bool            # Channel overlap flag
    segment_size::Union{Nothing,T}    # Segment size for averaging
    mean_rate1::Union{Nothing,T}      # Mean count rate of first signal
    mean_rate2::Union{Nothing,T}      # Mean count rate of second signal
end
```

## Constructor Functions

### Single-Segment Cross-Spectra

#### From Event Lists
```julia
CrossSpectrum(ev1::EventList, ev2::EventList, segment_size::Real, dt::Real;
              norm::String="frac", use_common_mean::Bool=true,
              fullspec::Bool=false, power_type::String="all") -> CrossSpectrum
```

**Parameters:**
- `ev1`, `ev2`: Input event lists containing photon arrival times
- `segment_size`: Duration of the analysis segment in seconds
- `dt`: Time bin size in seconds
- `norm`: Normalization scheme - `"frac"` (fractional RMS), `"leahy"` (Leahy), `"rms"` (RMS), `"abs"` (absolute), `"none"`
- `use_common_mean`: Use geometric mean of count rates for normalization
- `fullspec`: Include negative frequencies in output
- `power_type`: Power representation - `"all"` (complex), `"real"` (real part), `"absolute"` (magnitude)

**Usage:**


```julia
# Load NICER event data
event1 = readevents("ni1200120104_0mpu7_cl.evt", load_gti=true, sort=true)
ev_lowE  = filter_energy(x -> x >= 20 && x < 500, event1)
ev_highE = filter_energy(x -> x >= 500 && x <= 1500, event1)

# Basic cross-spectrum from 100-second segment with 10ms bins
cs = CrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01; norm="frac")
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 3199
      Normalization: frac
    




```julia
using Plots
plot(cs)
```

![](cross_spectrum_images\plot1.png)


```julia
# Leahy-normalized with negative frequencies
cs_leahy = CrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01, norm="leahy", fullspec=true)
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 6400
      Normalization: leahy
    



### Force linear scale when plotting full spectrum


```julia
plot(cs_leahy, xscale=:identity, yscale=:log10)
```

![](cross_spectrum_images\plot2.png)

#### from chandra data


```julia
# Load NICER event data
event1 = readevents("chandra_test.fits", load_gti=true, sort=true)
ev_lowE  = filter_energy(x -> x >= 20 && x < 500, event1)
ev_highE = filter_energy(x -> x >= 500 && x <= 1500, event1)

# Basic cross-spectrum from 100-second segment with 10ms bins
cs_chand = CrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01; norm="frac")
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 3199
      Normalization: frac
    




```julia
plot(cs_chand)
```

![](cross_spectrum_images\plot3.png)


```julia
cs_chand_leahy = CrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01, norm="leahy", fullspec=true)
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 6400
      Normalization: leahy
    




```julia
plot(cs_chand_leahy, xscale=:identity, yscale=:log10)
```

![](cross_spectrum_images\plot4.png)

#### From Light Curves
```julia
CrossSpectrum(lc1::LightCurve, lc2::LightCurve, segment_size::Union{Nothing,Real}=nothing;
              norm::String="frac", use_common_mean::Bool=true,
              fullspec::Bool=false, power_type::String="all") -> CrossSpectrum
```

**Parameters:**
- `lc1`, `lc2`: Pre-binned light curves
- `segment_size`: Analysis segment size (auto-determined if `nothing`)
- Other parameters identical to event list constructor

**Usage:**


```julia
event1 = readevents("ni1200120104_0mpu7_cl.evt", load_gti=true, sort=true)
ev_lowE  = filter_energy(x -> x >= 20 && x < 500, event1)
ev_highE = filter_energy(x -> x >= 500 && x <= 1500, event1)
# Create light curves
lc1 = create_lightcurve(ev_lowE, 0.01)
lc2 = create_lightcurve(ev_highE, 0.01)

# Instead of auto-determination, specify a reasonable segment size
# Use a segment size that's much smaller than your total observation time
total_time = maximum(lc1.time) - minimum(lc1.time)
segment_size = min(100.0, total_time / 10)  # 100 seconds or 1/10 of total time

println("Using segment size: ", segment_size, " seconds")
cs = CrossSpectrum(lc1, lc2, segment_size) ## or just go for "cs = CrossSpectrum(lc1, lc2, 200.0)"
```

    Using segment size: 100.0 seconds
    




    CrossSpectrum{Float64} (Single)
      Frequencies: 4999
      Normalization: frac
    




```julia
plot(cs)
```

![](cross_spectrum_images\plot5.png)


### Averaged Cross-Spectra

#### From Event Lists
```julia
AveragedCrossSpectrum(ev1::EventList, ev2::EventList, segment_size::Real, dt::Real;
                      norm::String="frac", use_common_mean::Bool=true,
                      fullspec::Bool=false, power_type::String="all",
                      fill_errors_on_creation::Bool=true) -> CrossSpectrum
```

**Usage:**


```julia
# Average multiple 64-second segments with 10ms binning
acs = AveragedCrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01)

# This creates a spectrum with m=53 segments averaged
println("Segments averaged: $(acs.m)")
println("Mean rates: $(acs.mean_rate1), $(acs.mean_rate2)")
```

    Segments averaged: 105
    Mean rates: 3078.7013392857143, 82.23645833333333
    


```julia
plot(acs)
```

![](cross_spectrum_images\plot6.png)

#### From Light Curves  
```julia
AveragedCrossSpectrum(lc1::LightCurve{T}, lc2::LightCurve{T}, segment_size::Real;
                      norm::String="frac", use_common_mean::Bool=true,
                      fullspec::Bool=false, power_type::String="all",
                      fill_errors_on_creation::Bool=true) -> CrossSpectrum
```

**Usage:**


```julia
# Average 100-second segments from light curves
acs = AveragedCrossSpectrum(lc1, lc2, 100.0)
```




    CrossSpectrum{Float64} (Averaged)
      Frequencies: 4999
      Segments averaged: 67
      Segment size: 100.0
      Mean rates: 3025.738208955224, 80.91417910447761
      Normalization: frac
    




```julia
plot(acs)
```

![](cross_spectrum_images\plot7.png)

## I guess u have noticed the difference between the plot while drawing with events and lightcurve, here is why

### 1. Data Processing Path
- **First plot (Light Curve method)**: Events → Light Curves → Cross Spectrum
- **Second plot (Event method)**: Events → Cross Spectrum directly

### 2. The Prominent Low-Frequency Features'

The "curves" or broad features you see in the first plot around 0.1-1 Hz are red noise characteristics that become more pronounced when using the light curve method. Here's why:

**Light Curve Method Effects:**
- Pre-binning the events into light curves can introduce correlations between adjacent bins
- GTI handling may create edge effects or gaps that manifest as low-frequency power
- The intermediate binning step can alter the noise characteristics, especially at low frequencies

**Event Method (Direct):**
- Processes raw photon arrival times directly
- More faithful to the original timing properties
- Better preserves the intrinsic noise characteristics

### 3. Physical Interpretation

The smoother, more curved appearance in the first plot suggests:
- **Stronger low-frequency power**: Indicating longer-timescale variability
- **Modified noise profile**: The pre-binning affects how noise propagates through the analysis
- **Potential artifacts**: Some of the low-frequency structure might be processing artifacts rather than astrophysical

## Utility Functions

### Spectrum Classification
```julia
is_averaged(cs::CrossSpectrum) -> Bool   # Returns true if m > 1
is_single(cs::CrossSpectrum) -> Bool     # Returns true if m = 1

println(is_averaged(cs))    # false (single segment)
println(is_averaged(acs))   # true (53 segments averaged)
```

### Noise Analysis

#### Theoretical Noise Level
```julia
theoretical_noise_level(cs::CrossSpectrum{T}) -> T
```
Calculate expected Poisson noise level based on count statistics.

**Usage:**


```julia
noise_level = theoretical_noise_level(acs)
println("Theoretical noise level: ", noise_level)
```

    Theoretical noise level: 0.0004938153701013313
    

#### White Noise Estimation

```julia
white_noise_level(cs::CrossSpectrum{T}; high_freq_fraction::Real=0.2) -> T
```
Estimate noise from high-frequency portion of spectrum (default: highest 20% of frequencies).

**Usage:**


```julia
empirical_noise = white_noise_level(cs, high_freq_fraction=0.3)
theoretical_noise = theoretical_noise_level(cs)
println("Empirical vs theoretical: ", empirical_noise, " vs ", theoretical_noise)
```

    Empirical vs theoretical: 0.0026649093667416293 vs 0.004056653203457984
    

#### Signal-to-Noise Ratio
```julia
signal_to_noise_ratio(cs::CrossSpectrum{T}) -> Vector{T}
```


```julia
snr = signal_to_noise_ratio(acs)
significant_mask = snr .> 3.0
println("Significant detections: ", sum(significant_mask))
```

    Significant detections: 886
    

#### to get direct all analysis use this


```julia
p1=plot(cs, Val(:noise_analysis))
p2=plot(acs, Val(:noise_analysis))
plot(p1, p2, layout=(1, 2), size=(1200, 500))
save_plot()
```

    Saved: powerspectrum_images/plot8.png
    




    "powerspectrum_images/plot8.png"



#### Error Estimation
```julia
fill_errors!(cs::CrossSpectrum{T}) -> CrossSpectrum{T}
```


```julia
# Populate power_err field with theoretical and sample variance
acs = AveragedCrossSpectrum(lc1, lc2, 100.0)
fill_errors!(acs)
# acs.power_err now contains error estimates
```




    CrossSpectrum{Float64} (Averaged)
      Frequencies: 4999
      Segments averaged: 67
      Segment size: 100.0
      Mean rates: 3025.738208955224, 80.91417910447761
      Normalization: frac
    



### Physical Properties

#### Coherence Function
```julia
coherence(cs::CrossSpectrum{T}) -> Vector{T}
```


```julia
# Calculate squared coherence γ²(f) = |Pxy|²/(Pxx·Pyy)
# Measures linear correlation vs frequency (0-1 scale)
coh = coherence(acs)
coh1 = coherence(cs)
p1 = plot(acs.freq, coh, xlabel="Frequency", ylabel="Coherence")
p2 = plot(acs.freq, coh1, xlabel="Frequency", ylabel="Coherence")
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot9.png)


```julia
p1=plot(acs, plot_type=:coherence) # or use `coh = coherence(acs) plot(acs.freq, coh, xlabel="Frequency", ylabel="Coherence", xscale=:log10)` 
p2=plot(cs, plot_type=:coherence)
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot10.png)

#### When using the plot method, it is on a log scale . you can either use this or  
```julia
coh = coherence(acs) 
plot(acs.freq, coh, xlabel="Frequency", ylabel="Coherence", xscale=:log10)

```

#### Phase and Time Lags
```julia
phase_lag(cs::CrossSpectrum) -> Vector    # Phase difference in radians
time_lag(cs::CrossSpectrum) -> Vector     # Time lag = phase/(2πf) in seconds
```


```julia
phase = phase_lag(cs)
tlag = time_lag(cs)
phase1 = phase_lag(acs)
tlag1 = time_lag(acs)
p1 = plot(cs.freq, phase, xlabel="Frequency (Hz)", ylabel="Phase Lag (rad)" )
p2 = plot(cs.freq, phase1, xlabel="Frequency (Hz)", ylabel="Phase Lag (rad)" )
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot11.png)


```julia
p1 = plot(cs, plot_type=:phase_lag, color=:red)
p2 = plot(acs, plot_type=:phase_lag, color=:red)
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot12.png)


```julia
p1 =plot(cs.freq, tlag, xlabel="Frequency (Hz)", ylabel="Time Lag (s)")
p2 =plot(cs.freq, tlag1, xlabel="Frequency (Hz)", ylabel="Time Lag (s)")
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot13.png)


```julia
p1 =plot(cs, plot_type=:time_lag, color=:orange)
p2 =plot(acs, plot_type=:time_lag, color=:orange)
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

![](cross_spectrum_images\plot14.png)

#### When using the plot method, it is on a log scale . you can either use this or  
```julia
plot(blah blah with, xscale=:log10)

```

## Rebinning Functions

### Linear Rebinning
```julia
rebin(cs::CrossSpectrum{T}, df_new::Real) -> CrossSpectrum{T}
rebin(cs::CrossSpectrum{T}, rebin_factor::Int) -> CrossSpectrum{T}


```julia
cs_rebinned = rebin(cs, 0.5)# Rebin to 0.5 Hz resolution
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 99
      Normalization: frac
    




```julia
cs_10x = rebin(cs, 10)          # Combine every 10 bins
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 499
      Normalization: frac
    




```julia
acs_rebinned = rebin(acs, 0.5)    # Rebin to 0.5 Hz resolution
```




    CrossSpectrum{Float64} (Averaged)
      Frequencies: 99
      Segments averaged: 67
      Segment size: 100.0
      Mean rates: 3025.738208955224, 80.91417910447761
      Normalization: frac
    




```julia
acs_10x = rebin(acs, 10)          # Combine every 10 bins
```




    CrossSpectrum{Float64} (Averaged)
      Frequencies: 499
      Segments averaged: 67
      Segment size: 100.0
      Mean rates: 3025.738208955224, 80.91417910447761
      Normalization: frac
    




```julia
p1 = plot(cs, cs_rebinned, Val(:rebinning_comparison), rebin_type="Linear")
p2 = plot(acs, acs_rebinned, Val(:rebinning_comparison), rebin_type="Linear")
plot(p1, p2, layout=(1, 2), size=(800, 400))
```

![](cross_spectrum_images\plot15.png)


```julia
p1 = plot(cs, cs_10x, Val(:rebinning_comparison), rebin_type="Linear")
p2 = plot(acs, acs_10x, Val(:rebinning_comparison), rebin_type="Linear")
plot(p1, p2, layout=(1, 2), size=(800, 400))
```

![](cross_spectrum_images\plot16.png)

### Logarithmic Rebinning
```julia
rebin_log(cs::CrossSpectrum{T}; f::Real=0.01) -> CrossSpectrum{T}


```julia
# Creates logarithmically spaced frequency bins
cs_log = rebin_log(cs, f=0.1)   # 10% fractional resolution
p1=plot(cs_log.freq, abs.(cs_log.power), xscale=:log10)
acs_log = rebin_log(acs, f=0.1)   # 10% fractional resolution
p2=plot(acs_log.freq, abs.(acs_log.power), xscale=:log10)
plot(p1, p2, layout=(1, 2), size=(1000, 400))
```

![](cross_spectrum_images\plot17.png)


```julia
cs_log_rebinned = rebin_log(cs, f=0.02)
p1=plot(cs, cs_log_rebinned, Val(:rebinning_comparison), rebin_type="Logarithmic")
acs_log_rebinned = rebin_log(acs, f=0.02)   
p2=plot(acs, acs_log_rebinned, Val(:rebinning_comparison), rebin_type="Logarithmic")
plot(p1, p2, layout=(1, 2), size=(1000, 400))
```

![](cross_spectrum_images\plot18.png)

### Geometric Rebinning
```julia
geometric_rebin(cs::CrossSpectrum{T}, factor::Real) -> CrossSpectrum{T}
```



```julia
# Each bin is 'factor' times wider than the previous
cs_geom_rebinned = geometric_rebin(cs, 1.3)
p1=plot(cs, cs_geom_rebinned, Val(:rebinning_comparison), rebin_type="Geometric")
acs_geom_rebinned = geometric_rebin(acs, 1.3)
p2=plot(acs, acs_geom_rebinned, Val(:rebinning_comparison), rebin_type="Geometric")
plot(p1, p2, layout=(1, 2), size=(1000, 400))
```

![](cross_spectrum_images\plot19.png)

### Adaptive Rebinning
```julia
adaptive_rebin(cs::CrossSpectrum{T}, target_snr::Real=3.0, 
               max_rebin_factor::Int=10) -> CrossSpectrum{T}


```julia
# Automatically rebin to achieve desired S/N ratio
cs_clean = adaptive_rebin(cs, 5.0)  # Rebin to achieve 5σ detection
```




    CrossSpectrum{Float64} (Single)
      Frequencies: 1249
      Normalization: frac
    




```julia
acs_clean = adaptive_rebin(acs, 5.0, 10)  # target_snr=5.0, max_rebin_factor=10
```




    CrossSpectrum{Float64} (Averaged)
      Frequencies: 2499
      Segments averaged: 67
      Segment size: 100.0
      Mean rates: 3025.738208955224, 80.91417910447761
      Normalization: frac
    




```julia
p1=plot(acs_clean)
p2=plot(cs_clean)
plot(p1, p2, layout=(1, 2), size=(1000, 400))
```

![](cross_spectrum_images\plot20.png)

## Diagnostic Functions

### Quality Assessment
```julia
quality_metrics(cs::CrossSpectrum{T}) -> Dict{String, Any}
noise_properties(cs::CrossSpectrum{T}) -> Dict{String, Any}


```julia
acs = AveragedCrossSpectrum(lc1, lc2, 100.0)
metrics = quality_metrics(acs)
println("Mean S/N: ", metrics["mean_snr"])
println("Significant fraction: ", metrics["significant_fraction"])

# If you want the actual number of significant detections:
total_bins = length(acs.freq)
significant_count = round(Int, metrics["significant_fraction"] * total_bins)
println("Significant detections: ", significant_count, " out of ", total_bins)
```

    Mean S/N: 19.328423543244085
    Significant fraction: 0.17723544708941788
    Significant detections: 886 out of 4999
    

### Significant Detection
```julia
significant_frequencies(cs::CrossSpectrum{T}, threshold::Real=3.0) -> Vector{T}
```


```julia
sig_freqs = significant_frequencies(acs, 5.0)  # Very conservative
println("Significant frequencies: ", length(sig_freqs), " out of ", length(acs.freq))
sig_freqs1 = significant_frequencies(cs, 5.0)  # Very conservative
println("Significant frequencies: ", length(sig_freqs1), " out of ", length(cs.freq))
```

    Significant frequencies: 630 out of 4999
    Significant frequencies: 354 out of 4999
    


```julia
p1=plot(acs, Val(:significant_detections), threshold=5.0)  # More conservative (5σ)
p2=plot(cs, Val(:significant_detections), threshold=5.0)
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot21.png)

### Aliasing Detection 
```julia
detect_aliasing(cs::CrossSpectrum{T}) -> (Bool, String)
```


```julia
aliased, message = detect_aliasing(cs)
if aliased
    println("Warning: ", message)
    # Consider using anti-aliasing filters or higher sampling rate
else
    println("not aliased")
end
```

    not aliased
    

## Visualization for more examples, please [visit here](https://github.com/StingraySoftware/Stingray.jl/pull/65#issue-3363831936)

### Basic Plotting 

Direct plotting with various types:

```julia
# Basic amplitude plot
plot(cs)
plot(acs)

# Or explicitly specify plot type
plot(acs, plot_type=:amplitude)
plot(cs, plot_type=:amplitude)
```

**Available Plot Types:**
- `:amplitude` - Cross-spectrum amplitude |Pxy(f)|
- `:power` - Cross-power |Pxy(f)|²  
- `:phase_lag` - Phase difference arg(Pxy(f))
- `:time_lag` - Time lag φ(f)/(2πf)
- `:coherence` - Squared coherence function
- `:snr` - Signal-to-noise ratio
- `:real_imaginary` - Real and imaginary components
- `:pds_comparison` - Input auto-power spectra

**Plot Customization:**
```julia
# With frequency range and color
plot(cs, plot_type=:amplitude, freq_range=(0.1, 10.0), color=:green)

# Phase lag with error bars
plot(acs, plot_type=:phase_lag, show_errors=true, freq_range=(0.1, 10))

# Coherence with noise level
plot(cs, plot_type=:coherence, show_noise_level=true)
```

### Advanced Visualization

#### Multi-Panel Analysis
```julia
# Comprehensive 2×3 analysis dashboard
plot(cs, Val(:analysis))
plot(acs, Val(:analysis))

# 2×2 noise diagnostic panel
plot(acs, Val(:noise_analysis))

# Aliasing and noise check
plot(cs, Val(:noise_diagnosis))
```

#### Specialized Plots
**Coherence with Confidence Levels:**


```julia
p1=plot(cs, Val(:coherence_confidence))
p2=plot(acs, Val(:coherence_confidence))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot22.png)

**Lag Analysis with Error Bars:**


```julia
# Time lag with reference values
p1=plot(cs, Val(:lag_frequency_errors), 
     freq_range=(1.0, 5.0),
     reference_lag=Dict(:frequency => 3.0, :expected_lag => 0.15))
p2=plot(acs, Val(:lag_frequency_errors), 
     freq_range=(1.0, 5.0),
     reference_lag=Dict(:frequency => 3.0, :expected_lag => 0.15))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot23.png)

### Phase frequency errors


```julia
p3=plot(cs, Val(:phase_frequency_errors),
     freq_range=(2.0, 5.0), 
     reference_phase=Dict(:frequency => 3.0, :expected_phase => π/3))
p4=plot(acs, Val(:phase_frequency_errors),
     freq_range=(2.0, 5.0), 
     reference_phase=Dict(:frequency => 3.0, :expected_phase => π/3))
plot(p3, p4, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot24.png)

### phase lag errors


```julia
p1=plot(cs, Val(:phase_lag_errors))
p2=plot(acs, Val(:phase_lag_errors))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot25.png)

### Time lag with rebinning


```julia
# Time lag with reference values
p1=plot(cs, Val(:lag_frequency_errors), rebin_factor=4, freq_range=(1.0, 10.0))
p2=plot(acs, Val(:lag_frequency_errors), rebin_factor=4, freq_range=(1.0, 10.0))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot26.png)

### **Signal-to-Noise Analysis:**


```julia
p1=plot(cs, Val(:frequency_snr), rebin_factor=5)
p2=plot(acs, Val(:noise_comparison))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot27.png)


```julia
p1=plot(cs, Val(:noise_timeline))
p2=plot(acs, Val(:noise_timeline))
plot(p1, p2, layout=(2, 1), size=(1000, 800))
```

![](cross_spectrum_images\plot28.png)

#### Comparison Plots

**Normalization Comparison:**
```julia
cs_leahy = CrossSpectrum(lc1, lc2, norm="leahy")
cs_frac = CrossSpectrum(lc1, lc2, norm="frac") 
cs_abs = CrossSpectrum(lc1, lc2, norm="abs")

plot([cs_leahy, cs_frac, cs_abs], Val(:normalization_comparison),
     norm_labels=["Leahy Norm", "Fractional RMS", "Absolute RMS"])
```

**Rebinning Effects:**
```julia
cs_rebinned = rebin(cs, 0.25)
plot(cs, cs_rebinned, Val(:rebinning_comparison), rebin_type="Linear")

cs_log_rebinned = rebin_log(acs, f=0.02)   
plot(acs, cs_log_rebinned, Val(:rebinning_comparison), rebin_type="Logarithmic")

cs_geom_rebinned = geometric_rebin(cs, 1.3)
plot(cs, cs_geom_rebinned, Val(:rebinning_comparison), rebin_type="Geometric")
```

## Normalization Schemes

### Fractional RMS Normalization (`norm="frac"`) - Default
Normalized by mean squared count rate. Measures fractional variability independent of source brightness. Best choice for most applications.

### Leahy Normalization (`norm="leahy"`)
Power spectrum normalized to have Poisson noise level of 2.0. Optimal for measuring absolute power levels and detecting periodic signals.

### Absolute RMS Normalization (`norm="rms"`)
Normalized to give RMS variability in count units. Preserves absolute amplitude information.

### Absolute Normalization (`norm="abs"`)
Minimal normalization preserving raw Fourier transform units.

### No Normalization (`norm="none"`)
Raw FFT output without correction factors.

## Best Practices

### Segment Size Selection
- **Single segments**: Use longest available continuous observation
- **Averaged spectra**: Balance frequency resolution vs statistical precision
- **Rule of thumb**: Segment size ≥ 10/f_min for reliable analysis at frequency f_min

### Normalization Choice
- **Fractional RMS**: General-purpose, source-brightness independent
- **Leahy**: Absolute power measurements, period detection
- **RMS**: When absolute variability amplitudes matter

### Error Analysis
- Always use `fill_errors_on_creation=true` for averaged spectra
- Check S/N ratios before interpreting results
- Use `quality_metrics()` to assess analysis reliability

### Rebinning Strategy
- **Linear**: Uniform frequency resolution needs
- **Logarithmic**: Power-law noise, broad-band features
- **Geometric**: Compromise between linear and logarithmic
- **Adaptive**: Automatic optimization for detection work

### Visualization Guidelines
- Use logarithmic scales for broad frequency ranges
- Show error bars for quantitative analysis  
- Include noise level overlays for context
- Apply frequency range limits to focus on features of interest

#### Normalization comparison (requires an array of differently normalized CrossSpectra)


```julia
cs_leahy = CrossSpectrum(lc1, lc2, norm="leahy")
cs_frac = CrossSpectrum(lc1, lc2, norm="frac") 
cs_abs = CrossSpectrum(lc1, lc2, norm="abs")

p1=plot([cs_leahy, cs_frac, cs_abs], Val(:normalization_comparison))

p2=plot([cs_leahy, cs_frac, cs_abs], Val(:normalization_comparison),
     norm_labels=["Leahy Norm", "Fractional RMS", "Absolute RMS"])

acs_leahy = AveragedCrossSpectrum(lc1, lc2,62, norm="leahy")
acs_frac = AveragedCrossSpectrum(lc1, lc2,62, norm="frac") 
acs_abs = AveragedCrossSpectrum(lc1, lc2,62, norm="abs")

p3=plot([acs_leahy, acs_frac, acs_abs], Val(:normalization_comparison))

p4=plot([acs_leahy, acs_frac, acs_abs], Val(:normalization_comparison),
     norm_labels=["Leahy Norm", "Fractional RMS", "Absolute RMS"])
plot(p1,p2,p3,p4, layout=(2, 2), size=(1200, 1200))
```

![](cross_spectrum_images\plot29.png)

## Example Workflow

```julia
# Load and filter event data
event1 = readevents("ni1050360108_0mpu7_cl.evt", load_gti=true, sort=true)
ev_lowE  = filter_energy(x -> x >= 20 && x < 500, event1)
ev_highE = filter_energy(x -> x >= 500 && x <= 1500, event1)

println("Low energy events: ", length(ev_lowE))
println("High energy events: ", length(ev_highE))

# Create light curves
lc1 = create_lightcurve(ev_lowE, 0.01)
lc2 = create_lightcurve(ev_highE, 0.01)

# Create both single and averaged cross-spectra
cs = CrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01; norm="frac")
acs = AveragedCrossSpectrum(ev_lowE, ev_highE, 64.0, 0.01)

# Quality assessment
metrics = quality_metrics(acs)
println("Mean S/N: $(metrics["mean_snr"])")
println("Significant detections: $(metrics["significant_fraction"]*100)%")

# Rebin for better S/N if needed
if metrics["mean_snr"] < 3.0
    acs_rebinned = adaptive_rebin(acs, target_snr=3.0)
else
    acs_rebinned = acs
end

# Comprehensive visualization
plot(acs_rebinned, Val(:analysis))

# Focus on coherence analysis
plot(acs_rebinned, Val(:coherence_confidence))

# Extract significant frequencies
sig_freqs = significant_frequencies(acs_rebinned, 3.0)
println("QPO candidates at: ", sig_freqs)

# Compare different rebinning approaches
cs_linear = rebin(acs, 0.25)
cs_log = rebin_log(acs, f=0.02)
cs_geom = geometric_rebin(acs, 1.3)

plot(acs, cs_linear, Val(:rebinning_comparison), rebin_type="Linear")
plot(acs, cs_log, Val(:rebinning_comparison), rebin_type="Logarithmic")
plot(acs, cs_geom, Val(:rebinning_comparison), rebin_type="Geometric")

# Noise analysis
plot(acs, Val(:noise_analysis))
```

## Troubleshooting

### Common Issues

#### 1. "No overlapping GTIs between event lists"
Check that both event lists cover the same time period and have valid GTI information.

#### 2. "No valid segments found for cross spectrum"  
Your segment size may be too large for the available data. Try smaller segments or check your GTI coverage.

#### 3. Poor Statistics/High Noise
Use adaptive rebinning or increase segment size for averaging to improve signal-to-noise ratio.

### Performance Tips

For large datasets:
- Use appropriate segment sizes (don't make them too small)
- Consider rebinning if you need many frequency bins
- Use `fill_errors_on_creation=false` initially if you don't need immediate error estimates

## thank you
