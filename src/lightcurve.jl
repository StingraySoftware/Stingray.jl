abstract type AbstractLightCurve{T} end

"""
    LightCurve{T}

A structure containing lightcurve data from a FITS file.

## Fields

- `timebins::Vector{T}`: Vector of time bins.
- `counts::Vector{Int}`: Vector of event counts in each time bin.
- `count_error::Vector{T}`: Vector of errors on the counts in each time bin.
- `err_method::Symbol`: Method used for computing the errors.
"""
struct LightCurve{T} <: AbstractLightCurve{T}
    timebins::Vector{T}
    counts::Vector{Int}
    count_error::Vector{T}
    err_method::Symbol
end

"""
   create_lightcurve(eventlist::EventList{T}, binsize::Real; err_method::Symbol=:poisson) -> LightCurve{T}

Create a lightcurve from an event list.

## Arguments
- `eventlist::EventList{T}`: Event list containing the event data.
- `binsize::Real`: Size of the time bins.

## Keyword Arguments
- `err_method::Symbol=:poisson`: Method for computing the errors.

## Returns
- [`LightCurve`](@ref) containing the binned data.
"""
function create_lightcurve(eventlist::EventList{T}, binsize::Real; err_method::Symbol=:poisson) where T
    # Validate error method first
    if err_method != :poisson
        throw(ArgumentError("Unsupported error computation method: $err_method"))
    end

    # Convert binsize to the same type as the event times
    binsize_t = convert(T, binsize)
    
    times = eventlist.times
    min_time = minimum(times)
    max_time = maximum(times)
    
    # Calculate number of bins needed
    time_range = max_time - min_time
    nbins = floor(Int, time_range / binsize_t)
    
    # If there's any remainder or if max_time falls exactly on a bin edge, add another bin
    if !isapprox(mod(time_range, binsize_t), zero(T)) || isapprox(mod(max_time - min_time, binsize_t), zero(T))
        nbins += 1
    end
    
    # Create arrays
    counts = zeros(Int, nbins)
    count_error = Vector{T}(undef, nbins)
    timebins = Vector{T}(undef, nbins + 1)
    
    # Create bin edges
    for i in 0:nbins
        timebins[i + 1] = min_time + (i * binsize_t)
    end
    
    # Bin the events using count() for better performance
    for i in 1:nbins
        bin_start = timebins[i]
        bin_end = timebins[i + 1]
        counts[i] = count(x -> bin_start <= x < bin_end, times)
    end
    
    # Calculate Poisson errors
    for i in 1:nbins
        count_error[i] = convert(T, sqrt(counts[i]))
    end
    
    return LightCurve{T}(collect(timebins), counts, count_error, err_method)
end

# Show method for pretty printing
function Base.show(io::IO, lc::LightCurve{T}) where T
    nbins = length(lc.counts)
    binsize = lc.timebins[2] - lc.timebins[1]
    print(io, "LightCurve{$T}(n=$nbins, binsize=$binsize)")
end