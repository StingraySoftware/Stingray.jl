```@meta
CurrentModule = Stingray
```

# Stingray

Documentation for [Stingray](https://github.com/matteobachetti/Stingray.jl).

A [Julia](http://julialang.org) package for X-ray Spectral Timing Analysis

This is a Julia porting of the Stingray software's 
[stingray](https://github.com/StingraySoftware/stingray) Python package.

## Usage Guide

Creating LightCurve from Time Stamps and Counts
```julia
julia> using Stingray

julia> times = range(0,49)

julia> counts = Int.(floor.(rand(50)*50000))

julia> lc = LightCurve(time=times,counts=counts, skip_checks=true, dt=1)
```

Plotting LightCurves
```julia
julia> using Plots

julia> plot(lc.time, lc.counts,xlabel="Time",ylabel="Counts")
```

Creating EventList from Photon Arrival Times

```julia
julia> times = [0.5, 1.1, 2.2, 3.7]

julia> mjdref=5800001

julia> ev = EventList(time=times, mjdref=mjdref)

julia> ev.energy = [0, 3, 4, 20.0]
```

Converting EventList to and from a LightCuvre

```julia
julia> ev = from_lc(lc)

julia> lc_new = to_lc(ev, 1)

julia> plot(lc.time, lc.counts-lc_new.counts)
```

Reading and Writing Events in a file
(Only FITS format supported currently)

Specify the filename and file format as parameters

```julia
julia> ev = EventList(time=time, mjdref=54000)

julia> write(ev, "file.fits", "fits")

julia> ev_new = read(EventList, "file.fits", "fits")
```
