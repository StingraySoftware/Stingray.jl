## Powerspectrum

Power spectrum analysis is the cornerstone of X-ray timing analysis, revealing the temporal variability characteristics of astrophysical sources. This guide demonstrates advanced usage of power spectrum functions across different X-ray missions, focusing on practical scientific applications.

### Key Concepts

| Concept | Definition | Scientific Impact |
|---------|------------|-------------------|
| **Power Spectrum** | Fourier transform of autocorrelation function | Reveals periodic signals, noise characteristics |
| **Normalization** | Mathematical scaling of power values | Enables comparison between different observations |
| **Averaging** | Combining multiple segments | Reduces statistical noise, improves S/N ratio |
| **Frequency Multiplication** | ν×P representation | Enhances high-frequency features, flattens red noise |

---

## Mission-Specific Considerations

Understanding each mission's capabilities is crucial for optimal power spectrum analysis:

| Mission | Optimal dt (s) | Typical Segment Size | Key Science |
|---------|----------------|---------------------|-------------|
| **Chandra** | 0.1 - 1.0 | 512 - 2048s | Long-term variability, eclipsing systems |
| **XMM** | 0.01 - 0.1 | 256 - 1024s | Intermediate timescales, energy-dependent timing |
| **RXTE** | 0.0001 - 0.001 | 64 - 256s | kHz QPOs, rapid variability |
| **NICER** | 0.001 - 0.01 | 128 - 512s | Neutron star oscillations, thermonuclear bursts |
| **NuSTAR** | 0.01 - 0.1 | 256 - 1024s | Hard X-ray timing, reflection lags |

---
## signatures:
```julia
# Single Power Spectrum from Light Curve
function Powerspectrum(lc::LightCurve{T}; norm::String = "frac") where {T<:Real}

# Averaged Power Spectrum from Light Curve
function AveragedPowerspectrum(
    lc::LightCurve{T},
    segment_size::Real;
    norm::String = "frac",
    epsilon::Real = 1e-5,
) where {T<:Real}

# Single Power Spectrum from Event List
function Powerspectrum(
    events::EventList{Vector{T},M},
    dt::Real,
    segment_size::Real;
    norm::String = "leahy",
) where {T<:Real,M}

# Averaged Power Spectrum from Event List
function AveragedPowerspectrum(
    events::EventList{Vector{T},M},
    segment_size::Real;
    norm::String = "frac",
    dt::Real = 1.0,
    epsilon::Real = 1e-5,
) where {T<:Real,M}
```
## 1. Data Loading and Mission Setup

Understanding your data characteristics before analysis is essential:


```julia
# Mission-optimized data loading
chandra = readevents("chandra_test.fits", mission="axaf", short=true)
xmm = readevents("xmm_test.fits", mission="xmm", short=true)
xte = readevents("xte_gx_test.evt.gz", mission="xte", short=true) 
nicer = readevents("ni1200120104_0mpu7_cl.evt", mission="nicer", short=true)
# nustar = readevents("nu90401330002A01_cl.evt", mission="nustar", short=true)

# Quick data assessment
println("Mission data loaded:")
for (name, events) in [("Chandra", chandra), ("XMM", xmm), ("RXTE", xte), ("NICER", nicer)]
    timespan = extrema(events.times)
    rate = length(events) / (timespan[2] - timespan[1])
    println("  $name: $(length(events)) events, $(round(rate, digits=1)) cts/s")
end
```

    Mission data loaded:
      Chandra: 4612 events, 4.9 cts/s
      XMM: 1708244 events, 1666.6 cts/s
      RXTE: 3518 events, 34.6 cts/s
      NICER: 21244574 events, 290.2 cts/s
    

**Why this matters**: Different missions have vastly different count rates and timing resolutions. RXTE can achieve μs timing, while Chandra excels at long-term monitoring.

---

## 2. Normalization Strategy

The choice of normalization fundamentally affects your scientific interpretation:

### Normalization Comparison Table

| Normalization | Formula | Units | Best For | Noise Level |
|---------------|---------|-------|----------|-------------|
| **Leahy** | 2N|ak|²/Ntot | dimensionless | QPO detection, χ² statistics | 2.0 |
| **Fractional** | 2dt|ak|²/μ²N | (rms/mean)² | Source comparison | 2dt/(μT) |
| **Absolute** | dt|ak|²/N | cts²·s | Physical power | dt·μ/T |


```julia
# normalization choices:
dt = 0.001  # 1ms resolution for neutron star timing
tstart, tstop = extrema(nicer.times)

# Create lightcurve with proper function signature
lc = create_lightcurve(nicer, dt, tstart=tstart, tstop=tstop)

# Leahy: Optimal for QPO detection and significance testing
ps_qpo = Powerspectrum(lc, norm="leahy")

# Fractional: Best for comparing different sources/observations
ps_compare = Powerspectrum(lc, norm="frac")

# abs: Physical interpretation of variability amplitude  
ps_amplitude = Powerspectrum(lc, norm="abs")

# Demonstrate the differences
println("Normalization effects at 100 Hz:")
freq_idx = argmin(abs.(ps_qpo.freq .- 100.0))
println("  Leahy: $(ps_qpo.power[freq_idx])")
println("  Fractional: $(ps_compare.power[freq_idx])") 
println("  abs: $(ps_amplitude.power[freq_idx])")
```

    Normalization effects at 100 Hz:
      Leahy: 0.3523331116572456
      Fractional: 0.001213911001828142
      abs: 102.26336311568573
    

# How abs and frac are related (useful sanity check)

Using the standard definitions:

- **Fractional rms normalization**:

$$
P_{\text{frac}} \propto \frac{2 \, \Delta t}{N \, \mu^2} \, |a_k|^2
$$

- **Absolute normalization**:

$$
P_{\text{abs}} \propto \frac{\Delta t}{N} \, |a_k|^2
$$

---

So,

$$
P_{\text{abs}} = P_{\text{frac}} \cdot \frac{\mu^2}{2}
$$

or equivalently,

$$
\mu = \sqrt{\frac{2 \, P_{\text{abs}}}{P_{\text{frac}}}}
$$

---

### Example with our numbers

Plugging in:

$$
\mu \approx \sqrt{\frac{2 \cdot 102.26336311568573}{0.001213911001828142}}
$$

$$
\mu \approx 410.5 \ \text{counts/s}
$$


### we have both methods you can either directly start analysis of powerspectrum and average powerspectrum via direct events we use [unbinned](https://github.com/StingraySoftware/Stingray.jl/blob/2c52daff3963bff886fad652c415d88c757e570d/src/gti.jl#L481) function there or via creating a light curve and then start analysis we use the [binned](https://github.com/StingraySoftware/Stingray.jl/blob/2c52daff3963bff886fad652c415d88c757e570d/src/gti.jl#L494) function for this :)

events[unbinned]


```julia
using Plots
```


```julia
events = readevents("ni1200120104_0mpu7_cl.evt", load_gti=true, sort=true)
ps = Powerspectrum(events, dt=10, norm="leahy")
```




    PowerSpectrum{Float64} (Single)
      Frequencies: 3659
      Frequency range: 1.4e-5 - 0.05 Hz
      Frequency resolution: 1.4e-5 Hz
      Power range: 0.278 - 3.7872147866e7
      Power errors: 2.0 - 2.0
      Normalization: leahy
      Total photons: 21208769
      Number of segments: 1
      Telescope: NICER
      Instrument: XTI
      Object: MAXI_J1820+070
      MJD reference: 0.0
      Time range: 1.3254001300020403e8 - 1.3261320799981469e8 s
      Bin size: 10.0 s
    




```julia
plot(ps)
```

![](powerspectrum_image\plot1.png)


```julia
avg_ps = AveragedPowerspectrum(events, 1024.0, norm="leahy", dt=0.1)
```




    AveragedPowerspectrum{Float64} (Averaged)
      Frequencies: 5119
      Frequency range: 0.000977 - 4.999 Hz
      Frequency resolution: 0.000977 Hz
      Segment size: 1024.0 s
      Segments averaged: 7
      Power range: 2.274 - 406329.574
      Power errors: 2.0 - 2.0
      Normalization: leahy
      Total photons: 19285536
      Mean rate: 2690.504 cts/s
      Telescope: NICER
      Instrument: XTI
      Object: MAXI_J1820+070
      MJD reference: 0.0
      Time range: 1.3254001300020403e8 - 1.3261320799981469e8 s
      Bin size: 0.1 s
      GTI intervals: 16
      Total good time: 8285.455 s
      Analysis method: direct_events_processing
      Original events: 21244574
      Time resolution: 0.1 s
    




```julia
plot(avg_ps)
```

![](powerspectrum_image\plot2.png)

lightcurve[binned]


```julia
lc = create_lightcurve(events, 10)
ps_lc = Powerspectrum(lc, norm="leahy")
```




    PowerSpectrum{Float64} (Single)
      Frequencies: 3659
      Frequency range: 1.4e-5 - 0.05 Hz
      Frequency resolution: 1.4e-5 Hz
      Power range: 1.313 - 3.7918349935e7
      Power errors: 2.0 - 2.0
      Normalization: leahy
      Total photons: 21244574
      Number of segments: 1
      Telescope: NICER
      Instrument: XTI
      Object: MAXI_J1820+070
      MJD reference: 56658.0
      Time range: 1.3254001300020403e8 - 1.3261320799981469e8 s
      Bin size: 10.0 s
    




```julia
plot(ps_lc)
```

![](powerspectrum_image\plot3.png)


```julia
ps_avg_lc = AveragedPowerspectrum(lc, 1024.0, norm="leahy")
```




    AveragedPowerspectrum{Float64} (Averaged)
      Frequencies: 50
      Frequency range: 0.00098 - 0.049 Hz
      Frequency resolution: 0.00098 Hz
      Segment size: 1024.0 s
      Segments averaged: 7
      Power range: 882.054 - 364567.336
      Power errors: 2.0 - 2.0
      Normalization: leahy
      Total photons: 19544935
      Mean rate: 2726.693 cts/s
      Telescope: NICER
      Instrument: XTI
      Object: MAXI_J1820+070
      MJD reference: 56658.0
      Time range: 1.3254001300020403e8 - 1.3261320799981469e8 s
      Bin size: 10.0 s
      GTI intervals: 16
      Total good time: 8285.455 s
    




```julia
plot(ps_avg_lc)
```

![](powerspectrum_image\plot4.png)

## 3. Averaging Strategies

Segment size selection is a critical trade-off between frequency resolution and statistical precision:

### Segment Size Trade-offs

| Segment Size | Frequency Resolution | Statistical Error | Best Science Case |
|--------------|---------------------|-------------------|-------------------|
| **Short (64s)** | Δf = 0.016 Hz | High noise | Broad-band noise, kHz QPOs |
| **Medium (256s)** | Δf = 0.004 Hz | Moderate | mHz QPOs, intermediate timescales |
| **Long (1024s)** | Δf = 0.001 Hz | Low noise | μHz phenomena, long-term trends |


```julia
# kHz QPO search: Short segments, high time resolution
ps_khz = AveragedPowerspectrum(events, 64.0, norm="leahy", dt=0.001)
println("kHz QPO setup: Δf=$(ps_khz.df) Hz, fmax=$(maximum(ps_khz.freq)) Hz")

# Low-frequency QPO: Long segments for frequency precision  
ps_lfqpo = AveragedPowerspectrum(events, 1024.0, norm="leahy", dt=0.001)
println("LFQPO setup: Δf=$(ps_lfqpo.df) Hz, $(ps_lfqpo.m) segments")

# Broad-band noise characterization: Balanced approach
ps_broadband = AveragedPowerspectrum(events, 256.0, norm="frac", dt=0.001)
println("Broadband: $(ps_broadband.m) segments, σ_rel=$(2.0/sqrt(ps_broadband.m))")
```

    kHz QPO setup: Δf=0.015625 Hz, fmax=499.984375 Hz
    LFQPO setup: Δf=0.0009765625 Hz, 7 segments
    Broadband: 30 segments, σ_rel=0.3651483716701107
    

### Do check for Total duration


```julia
t_start = minimum(events.times)
t_end   = maximum(events.times)
total_T = t_end - t_start
println("Total duration = $total_T seconds")
```

    Total duration = 73194.99961066246 seconds
    

**Critical insight**: The number of segments (m) determines your statistical precision. Error scales as 1/√m, so 4× more segments = 2× better precision.

---

## 4. Plotting and Visualization 

to see all possible plots may consider visiting [here](https://github.com/StingraySoftware/Stingray.jl/pull/59#issue-3273135469)

## Basic Plotting Functions

### Single Power Spectrum Plot
```julia
using Plots

# Basic plot
plot(ps_avg)

# With noise level
plot(ps_avg, show_noise=true)

# With error bars  
plot(ps_avg, show_errors=true, show_noise=true)

# Frequency-multiplied (ν×P vs P) - enhances high-frequency features
plot(ps_avg, freq_mult=true, show_noise=true)

# Linear scale (rare, usually log is better)
plot(ps_avg, log_scale=false)

# Custom axis limits
plot(ps_avg, axis_limits=[0.01, 100, 1e-5, 1e-1])  # [fmin, fmax, pmin, pmax]
plot(ps_avg, axis_limits=[0.1, 50])                # [fmin, fmax] only
```
### QPO Detection Protocol

QPO detection requires specific plotting strategies and significance assessment:


```julia
ps = AveragedPowerspectrum(xte, 8.0, norm="leahy", dt=0.0001)
# Standard QPO search plot with significance levels
plot(ps, Val(:qpo),
     qpo_range=(0.1, 2000),         # Cover both LFQPO and kHz QPO ranges
     freq_mult=true,                # Flatten red noise continuum  
     significance_level=3.0,        # 3σ detection threshold
     mark_peaks=true)
```

![](powerspectrum_image\plot5.png)


```julia
# Multi-scale QPO analysis
p1 = plot(ps, Val(:qpo), qpo_range=(0.01, 10), title="Low-Frequency QPOs")
```

![](powerspectrum_image\plot6.png)


```julia
p2 = plot(ps, Val(:qpo), qpo_range=(10, 100), title="Intermediate QPOs") 
```

![](powerspectrum_image\plot7.png)


```julia
p3 = plot(ps, Val(:qpo), qpo_range=(100, 2000), title="kHz QPOs")
```

![](powerspectrum_image\plot8.png)

### Noise Subtraction for Intrinsic Variability and using `freq_mult=true`


```julia
# Remove Poisson noise to study intrinsic source variability
ps_intrinsic = AveragedPowerspectrum(nicer, 256.0, norm="frac", dt=0.01)

plot(ps_intrinsic, 
     subtract_noise=true,           # Remove Poisson component
     show_noise=true,               # Show original noise level for reference
     freq_mult=true,
     title="Intrinsic Variability (Poisson Subtracted)")
```

![](powerspectrum_image\plot9.png)

---

## 5. Multi-Band Analysis

### Energy-Dependent Timing Analysis

Energy-resolved timing reveals spectral-timing correlations crucial for understanding emission mechanisms:


```julia
lc = create_lightcurve(events, 0.1)
events_soft = filter_energy(e -> 30 <= e <= 200, events)    # ~0.3-2 keV
events_medium = filter_energy(e -> 200 <= e <= 800, events) # ~2-8 keV  
events_hard = filter_energy(e -> 800 <= e <= 1200, events)  # ~8-12 keV

# Now this should work:
ps_soft = AveragedPowerspectrum(events_soft, 256.0, norm="frac", dt=0.1)
ps_medium = AveragedPowerspectrum(events_medium, 256.0, norm="frac", dt=0.1)
ps_hard = AveragedPowerspectrum(events_hard, 256.0, norm="frac", dt=0.1)

plot([ps_soft, ps_medium, ps_hard], 
     labels=["Soft", "Medium", "Hard"], 
     colors=[:blue, :green, :red])
```

![](powerspectrum_image\plot10.png)

using `freq_mult=true`


```julia
plot([ps_soft, ps_medium, ps_hard], segment_size=256.0, norm="frac",
     colors=[:blue, :green, :red], freq_mult=true)
```

![](powerspectrum_image\plot11.png)

---

## 6. Advanced Statistical Analysis

### Error Analysis and Confidence Intervals


```julia
# Comprehensive error analysis
ps = AveragedPowerspectrum(nicer, 256.0, norm="leahy", dt=0.1)

# Extract statistical properties
frequencies = freqs(ps)
power_values = power(ps) 
power_errors = errors(ps)

# Calculate confidence intervals (for Leahy normalization)
conf_95_lower = power_values .- 1.96 * power_errors
conf_95_upper = power_values .+ 1.96 * power_errors

println("Statistical Summary:")
println("  Segments averaged: $(ps.m)")
println("  Expected noise level: 2.0 (Leahy)")
println("  Observed high-freq mean: $(mean(power_values[end-50:end]))")
println("  Theoretical uncertainty: $(2.0/sqrt(ps.m))")

# Plot with confidence intervals
plot(ps, show_errors=true, error_alpha=0.3,
     title="Power Spectrum with 95% Confidence Intervals")
```

    Statistical Summary:
      Segments averaged: 30
      Expected noise level: 2.0 (Leahy)
      Observed high-freq mean: 11.00980554096393
      Theoretical uncertainty: 0.3651483716701107
    

![](powerspectrum_image\plot12.png)

### Dead Time Corrections


```julia
# Mission-specific dead time corrections
function apply_deadtime_correction(ps, mission)
    if mission == "xte"
        dead_time = 10e-6  # 10 μs for RXTE/PCA
        correction = 1.0 ./ (1.0 .- 2.0 * dead_time * ps.freq)
        return ps.power .* correction
    elseif mission == "nicer"  
        dead_time = 2.5e-6  # 2.5 μs for NICER
        correction = 1.0 ./ (1.0 .- 2.0 * dead_time * ps.freq)
        return ps.power .* correction
    end
    return ps.power  # No correction needed
end

# Calculate correction factor separately for visualization
function get_correction_factor(ps, mission)
    if mission == "xte"
        dead_time = 10e-6
        return 1.0 ./ (1.0 .- 2.0 * dead_time * ps.freq)
    elseif mission == "nicer"  
        dead_time = 2.5e-6
        return 1.0 ./ (1.0 .- 2.0 * dead_time * ps.freq)
    end
    return ones(length(ps.freq))
end

# Create corrected powerspectrum
corrected_power = apply_deadtime_correction(ps, "nicer")
ps_corrected = AveragedPowerspectrum(
    ps.freq, corrected_power, ps.power_err, ps.norm,
    ps.df, ps.segment_size, ps.nphots, ps.m, 
    ps.mean_rate, ps.n, ps.metadata
)
correction_factor = get_correction_factor(ps, "nicer")
fractional_diff = (corrected_power - ps.power) ./ ps.power * 100  # in percent
high_freq_mask = ps.freq .> 1000  # Only frequencies > 1 kHz
if sum(high_freq_mask) > 0
    plot([ps.power[high_freq_mask], corrected_power[high_freq_mask]], 
         ps.freq[high_freq_mask],
         labels=["Original", "Corrected"],
         xscale=:log10,
         yscale=:log10,
         xlabel="Frequency (Hz)",
         ylabel="Power",
         title="Dead Time Correction (High Frequencies Only)")
end
ratio = corrected_power ./ ps.power
println("Dead time correction statistics:")
println("Frequency range: $(minimum(ps.freq)) - $(maximum(ps.freq)) Hz")
println("Correction factor range: $(minimum(correction_factor)) - $(maximum(correction_factor))")
println("Maximum fractional change: $(maximum(fractional_diff))%")
println("At 1 kHz: $(correction_factor[findmin(abs.(ps.freq .- 1000))[2]]))")
println("At 5 kHz: $(correction_factor[findmin(abs.(ps.freq .- 5000))[2]]))")
```

    Dead time correction statistics:
    Frequency range: 0.00390625 - 4.99609375 Hz
    Correction factor range: 1.0000000195312504 - 1.0000249810927895
    Maximum fractional change: 0.002498109278943023%
    At 1 kHz: 1.0000249810927895)
    At 5 kHz: 1.0000249810927895)
    


```julia
plot(ps.freq, correction_factor, 
     xscale=:log10, 
     xlabel="Frequency (Hz)", 
     ylabel="Correction Factor",
     title="NICER Dead Time Correction Factor",
     label="Correction = 1/(1-2τf)")
```

![](powerspectrum_image\plot13.png)

---

## Key Performance Guidelines

### Computational Efficiency

| Data Size | Recommended Approach | Memory Usage | Processing Time |
|-----------|---------------------|--------------|-----------------|
| < 10⁶ events | Direct processing | ~100 MB | < 1 minute |
| 10⁶ - 10⁷ events | Chunked analysis | ~500 MB | 1-10 minutes |
| > 10⁷ events | Parallel processing | ~2 GB | 10+ minutes |

# Complete X-ray Power Spectrum Analysis Guide

## Overview
This guide demonstrates all power spectrum functions using real X-ray mission data. Power spectra are fundamental tools in X-ray astronomy for studying variability, detecting QPOs (Quasi-Periodic Oscillations), and characterizing noise properties.

---

## 1. Basic Data Loading (EventList Creation)

### Loading Different Mission Data
```julia
# Load data from different X-ray missions
chandra = readevents("chandra_test.fits", mission="axaf", short=true)
xmm = readevents("xmm_test.fits", mission="xmm", short=true)  
xte = readevents("xte_gx_test.evt.gz", mission="xte", short=true)
nicer = readevents("ni1200120104_0mpu7_cl.evt", mission="nicer", short=true)
nustar = readevents("nustar_fpma_cl.evt", mission="nustar", short=true)
swift = readevents("swift_xrt.evt", mission="swift", short=true)

# Basic EventList information
println("Chandra: $(length(chandra)) events")
println("XMM: $(length(xmm)) events") 
println("RXTE: $(length(xte)) events")
println("NICER: $(length(nicer)) events")
```

---

## 2. Single Power Spectrum (PowerSpectrum)

### Purpose: Analyze variability in a single time segment

```julia
# From Light Curve
lc = LightCurve(chandra, dt=0.1)  # 0.1s bins
ps_single = Powerspectrum(lc, norm="leahy")

# From EventList directly  
ps_events = Powerspectrum(chandra, dt=0.1, segment_size=1024.0, norm="leahy")

# Different normalizations - each has specific uses:
ps_leahy = Powerspectrum(lc, norm="leahy")      # For Poisson statistics
ps_frac = Powerspectrum(lc, norm="frac")        # Fractional variability  
ps_rms = Powerspectrum(lc, norm="rms")          # RMS variability
ps_abs = Powerspectrum(lc, norm="abs")          # Absolute power

# Access properties
println("Frequencies: ", ps_single.freq[1:5])   # First 5 frequencies
println("Power values: ", ps_single.power[1:5]) # First 5 power values
println("Normalization: ", ps_single.norm)      # Normalization type
println("Frequency resolution: ", ps_single.df) # Δf
```

---

## 3. Averaged Power Spectrum (AveragedPowerspectrum)

### Purpose: Reduce noise by averaging multiple segments

```julia
# From Light Curve
lc_long = LightCurve(nicer, dt=0.01)  # High time resolution
ps_avg = AveragedPowerspectrum(lc_long, 256.0, norm="frac")

# From EventList (recommended for precision)
ps_avg_events = AveragedPowerspectrum(nicer, 256.0, norm="frac", dt=0.01)

# Different segment sizes for different science goals:
ps_short = AveragedPowerspectrum(nicer, 64.0, norm="frac", dt=0.01)   # High freq resolution
ps_medium = AveragedPowerspectrum(nicer, 256.0, norm="frac", dt=0.01) # Balanced
ps_long = AveragedPowerspectrum(nicer, 1024.0, norm="frac", dt=0.01)  # Low noise

# Properties specific to averaged spectra
println("Number of segments: ", ps_avg.m)
println("Mean count rate: ", ps_avg.mean_rate, " cts/s")
println("Segment size: ", ps_avg.segment_size, " s")
```

---



### EventList Direct Plotting
```julia
# Plot directly from events (creates power spectrum internally)
plot(nicer, segment_size=256.0, norm="frac", show_noise=true)

# High time resolution for timing analysis
plot(nicer, segment_size=128.0, dt=0.001, norm="leahy", freq_mult=true)

# Different normalizations
plot(xte, segment_size=512.0, norm="rms", show_errors=true)
```

---

## 5. Comparing Multiple Power Spectra

### Vector of Power Spectra
```julia
# Create multiple spectra with different segment sizes
ps1 = AveragedPowerspectrum(nicer, 128.0, norm="frac", dt=0.01)
ps2 = AveragedPowerspectrum(nicer, 256.0, norm="frac", dt=0.01) 
ps3 = AveragedPowerspectrum(nicer, 512.0, norm="frac", dt=0.01)

spectra_vector = [ps1, ps2, ps3]

# Plot comparison
plot(spectra_vector, 
     labels=["128s segments", "256s segments", "512s segments"],
     colors=[:blue, :red, :green],
     show_noise=true,
     freq_mult=true)

# Different missions comparison
ps_chandra = AveragedPowerspectrum(chandra, 256.0, norm="frac", dt=0.1)
ps_xmm = AveragedPowerspectrum(xmm, 256.0, norm="frac", dt=0.1)
ps_xte = AveragedPowerspectrum(xte, 256.0, norm="frac", dt=0.001)

mission_spectra = [ps_chandra, ps_xmm, ps_xte]
plot(mission_spectra,
     labels=["Chandra", "XMM", "RXTE"],
     alpha=0.8,
     linewidth=2)
```

---

## 6. Multi-Band Spectral-Timing Analysis

### Purpose: Study energy-dependent variability

```julia
# Filter events by energy bands
nicer_soft = filter_energy(nicer, 0.5, 2.0)    # Soft X-rays
nicer_hard = filter_energy(nicer, 2.0, 10.0)   # Hard X-rays

# Create dictionary for multi-band analysis
events_dict = Dict(
    "0.5-2 keV" => nicer_soft,
    "2-10 keV" => nicer_hard
)

# Plot multi-band comparison
plot(events_dict, 
     segment_size=256.0, 
     energy_labels=["Soft", "Hard"],
     colors=[:blue, :red],
     freq_mult=true)

# More energy bands (if mission supports it)
if mission == "nustar" || mission == "xmm"
    events_multiband = Dict(
        "3-5 keV" => filter_energy(events, 3.0, 5.0),
        "5-10 keV" => filter_energy(events, 5.0, 10.0), 
        "10-20 keV" => filter_energy(events, 10.0, 20.0),
        "20-40 keV" => filter_energy(events, 20.0, 40.0)
    )
    
    plot(events_multiband,
         segment_size=512.0,
         colors=[:blue, :green, :orange, :red])
end
```

---

## 7. QPO Analysis (Specialized Plotting)

### Purpose: Search for Quasi-Periodic Oscillations

```julia
# QPO-focused analysis
ps_qpo = AveragedPowerspectrum(xte, 256.0, norm="leahy", dt=0.001)

# QPO search plot
plot(ps_qpo, Val(:qpo),
     qpo_range=(0.1, 100),           # Frequency range to search
     freq_mult=true,                 # ν×P enhances QPO detection
     significance_level=3.0,         # 3σ detection threshold
     mark_peaks=true)                # Mark potential QPOs

# Different QPO frequency ranges
plot(ps_qpo, Val(:qpo), qpo_range=(0.01, 10))   # Low-frequency QPOs  
plot(ps_qpo, Val(:qpo), qpo_range=(10, 1000))   # High-frequency QPOs

# With custom noise level
mean_rate = mean_countrate(xte)
noise_level = poisson_level("leahy", meanrate=mean_rate)
plot(ps_qpo, Val(:qpo), 
     noise_level=noise_level,
     significance_level=5.0)         # 5σ for conservative detection
```

---


---

## 9. Advanced Analysis Workflows

### Scientific Use Cases

```julia
# 1. TIMING ANALYSIS WORKFLOW
# ==========================
events = readevents("accreting_bh.evt", mission="nicer", short=true)

# High time resolution for rapid variability
ps_fast = AveragedPowerspectrum(events, 64.0, norm="rms", dt=0.0001)
plot(ps_fast, freq_mult=true, show_noise=true, axis_limits=[1, 4000])

# 2. QPO DETECTION WORKFLOW  
# ==========================
# Multiple segment sizes for different QPO types
ps_lfqpo = AveragedPowerspectrum(events, 512.0, norm="leahy", dt=0.001)  # Low-freq QPOs
ps_hfqpo = AveragedPowerspectrum(events, 128.0, norm="leahy", dt=0.0001) # High-freq QPOs

plot(ps_lfqpo, Val(:qpo), qpo_range=(0.01, 30))
plot(ps_hfqpo, Val(:qpo), qpo_range=(30, 2000))

# 3. ENERGY-DEPENDENT TIMING
# ===========================
soft_events = filter_energy(events, 0.5, 2.0)
hard_events = filter_energy(events, 2.0, 10.0)

ps_soft = AveragedPowerspectrum(soft_events, 256.0, norm="frac")
ps_hard = AveragedPowerspectrum(hard_events, 256.0, norm="frac")

# Compare lag properties
plot([ps_soft, ps_hard], 
     labels=["Soft", "Hard"], 
     freq_mult=true)

# 4. NOISE CHARACTERIZATION
# =========================
# Subtract Poisson noise to study intrinsic variability
plot(ps, subtract_noise=true, show_noise=true, freq_mult=true)

# 5. MISSION COMPARISON
# ====================
ps_nicer = AveragedPowerspectrum(nicer, 256.0, norm="frac", dt=0.01)
ps_xte = AveragedPowerspectrum(xte, 256.0, norm="frac", dt=0.001)

plot([ps_nicer, ps_xte],
     labels=["NICER", "RXTE"],
     alpha=0.7,
     axis_limits=[0.01, 100, 1e-6, 1e-1])
```

---

## 10. Error Analysis and Statistics

```julia
ps = AveragedPowerspectrum(events, 256.0, norm="leahy")

# Statistical properties
println("Number of segments averaged: $(ps.m)")
println("Expected Leahy noise level: 2.0")
println("Actual mean power at high freq: $(mean(ps.power[end-10:end]))")

# Error propagation for different normalizations
if ps.norm == "leahy"
    println("Leahy errors: χ² distributed, σ = 2")
elseif ps.norm in ["frac", "rms"]  
    println("Fractional errors scale as 1/√m where m = $(ps.m)")
end

# Plot with error bars for statistical assessment
plot(ps, show_errors=true, error_alpha=0.5, show_noise=true)
```

---

## 11. Custom Analysis Parameters

```julia
# Fine-tuning for specific science goals

# ULTRA-HIGH TIME RESOLUTION (for kHz QPOs)
ps_ultra = AveragedPowerspectrum(events, 32.0, norm="leahy", dt=0.0001)
plot(ps_ultra, axis_limits=[100, 4000], freq_mult=true)

# LONG-TERM VARIABILITY (low frequencies)
ps_longterm = AveragedPowerspectrum(events, 2048.0, norm="frac", dt=1.0)
plot(ps_longterm, axis_limits=[1e-4, 1], show_noise=true)

# CUSTOM PLOTTING STYLES
plot(ps, 
     drawstyle=:steppost,        # Step function (default)
     linewidth=2.0,
     noise_style=:dot,           # Dotted noise line
     error_alpha=0.3)            # Semi-transparent errors

# PUBLICATION-READY PLOTS
plot(ps,
     show_noise=true,
     show_errors=true, 
     freq_mult=true,
     axis_limits=[0.01, 100, 1e-5, 1e-1],
     title="NICER Power Spectrum",
     xlabel="Frequency (Hz)",
     ylabel="Fractional Power × ν")
```

---

## Summary of Function Purposes

| Function | Purpose | When to Use |
|----------|---------|-------------|
| `Powerspectrum()` | Single segment analysis | Quick look, short data |
| `AveragedPowerspectrum()` | Multi-segment averaging | Reduce noise, long observations |
| `plot(ps)` | Basic visualization | Initial data inspection |
| `plot(ps, show_noise=true)` | Compare with Poisson level | Assess significance |
| `plot(ps, freq_mult=true)` | Enhance high frequencies | QPO detection |
| `plot(vector_of_ps)` | Compare multiple spectra | Parameter studies |
| `plot(events_dict)` | Energy-dependent timing | Spectral-timing analysis |
| `plot(ps, Val(:qpo))` | QPO search | Periodic signal detection |
| `freqs()`, `power()`, `errors()` | Data extraction | Custom analysis |

### Key Normalizations:
- **"leahy"**: Poisson statistics, mean=2 for pure noise
- **"frac"**: Fractional RMS variability  
- **"rms"**: RMS variability in same units as light curve
- **"abs"**: Absolute power (rarely used)

This guide covers all functions with practical X-ray astronomy applications!

## Thank you
