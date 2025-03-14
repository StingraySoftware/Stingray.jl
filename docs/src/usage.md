```@meta
CurrentModule = Stingray
```

# Stingray

Documentation for [Stingray](https://github.com/matteobachetti/Stingray.jl).

```@index
```

```@autodocs
Modules = [Stingray]
```

## Installing Julia
Before installing Stingray.jl, ensure you have **Julia v1.7 or later** installed. If you haven’t installed Julia yet, follow these steps:
1. Download Julia from the official site: [https://julialang.org/downloads/](https://julialang.org/downloads/)
2. Follow the installation instructions for your operating system.
3. Add Julia to your system’s PATH for easy command-line access.

## Installing Stingray.jl
To install Stingray.jl, open a Julia REPL (press `]` to enter the package manager) and run:

```julia
using Pkg
Pkg.add("Stingray")
```

This will download and install all necessary dependencies.

## Cloning the Development Version
If you want the latest features, clone the GitHub repository and work in development mode:

```julia
Pkg.add(url="https://github.com/matteobachetti/Stingray.jl")
Pkg.develop("Stingray")
```

## Checking Your Installation
To verify the installation, run the following in the Julia REPL:

```julia
using Stingray
@info "Stingray.jl loaded successfully!"
```

# Usage Guide for Stingray.jl

Documentation for [Stingray.jl](https://github.com/matteobachetti/Stingray.jl).

```@index
```

## Overview
Stingray.jl provides powerful tools for **X-ray spectral timing analysis**, including **Fourier transforms, periodograms, and coherence calculations**. This guide covers the basic usage of the package.

## Importing the Package
After installation, load Stingray.jl in your Julia session:
```julia
using Stingray
```

## Creating a Simulated Light Curve
```julia
using Random, Plots

# Generate Simulated Light Curve
N = 1024  # Number of data points
t = collect(0:0.1:(N-1)*0.1)  # Time array
Random.seed!(42)
light_curve = sin.(2π * 0.5 .* t) + 0.3 * randn(N)  # Sine wave + noise

# Plot Light Curve
plot(t, light_curve, label="Simulated Light Curve", xlabel="Time (s)", ylabel="Intensity", title="X-ray Light Curve", legend=:topright)
```

## Computing a Power Spectrum
```julia
using FFTW

# Compute Fourier Transform
frequencies = fftshift(fftfreq(N, 0.1))  # Compute frequency bins
power_spectrum = abs.(fftshift(fft(light_curve))).^2  # Compute power spectrum

# Plot Power Spectrum
plot(frequencies, power_spectrum, xlabel="Frequency (Hz)", ylabel="Power", title="Power Spectrum", legend=false)
```

```@autodocs
Modules = [Stingray]
```

## Power Colors

### Overview

The `power_color` function is a key component in Stingray.jl, designed to analyze frequency-power distributions.

### Function: `power_color`

```@autodocs
Modules = [Stingray]
Pages = ["power_color.jl"]
```

### **Description**
The `power_color` function calculates power color indices based on frequency and power data.

### **Usage**
```julia
pc0, pc0e, pc1, pc1e = power_color(freq, power)
```

### **Arguments**
- `freq::Vector{Float64}`: The frequency values.
- `power::Vector{Float64}`: The corresponding power values.
- `return_log::Bool=false` (optional): Returns logarithmic power color if `true`.

### **Returns**
- `pc0::Float64`: Power color index 0.
- `pc0e::Float64`: Error in pc0.
- `pc1::Float64`: Power color index 1.
- `pc1e::Float64`: Error in pc1.

### **Example**
```julia
freq = collect(0.0001:0.00001:17)
power = 1 ./ freq
pc0, pc0e, pc1, pc1e = power_color(freq, power)
println("Power Color Indices: ", pc0, pc0e, pc1, pc1e)
```

# Temporary Reference
# Refer to this repository for more content: [GitHub - Notebooks](https://github.com/kashish2210/notebooks)
