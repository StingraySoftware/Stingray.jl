"""
    LightCurve{T}

A structure representing a light curve, which is a time series of event counts.

## Fields

- timebins::Vector{T}: Time bins for the light curve.
- counts::Vector{Int}: Number of events in each time bin.
- count_error::Vector{T}: Error estimate for each bin.
- err_method::Symbol: Method used for error estimation (:poisson).
"""
struct LightCurve{T}
    timebins::Vector{T}
    counts::Vector{Int}
    count_error::Vector{T}
    err_method::Symbol
    
    function LightCurve{T}(timebins, counts, count_error, err_method) where T
        length(counts) == length(count_error) || throw(ArgumentError("Counts and error arrays must have the same length."))
        return new(timebins, counts, count_error, err_method)
    end
end

"""
    create_lightcurve(eventlist::EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T

Generate a light curve from an EventList by binning event times.

## Arguments

- eventlist::EventList{T}: The input event list.
- bin_size::T: The size of each time bin.
- err_method::Symbol: The method for error estimation (:poisson).

## Returns

A LightCurve instance containing binned event counts.

## Notes

- The function bins event times into intervals of bin_size.
- Errors are estimated using Poisson statistics.
- Handles cases where no events fall within a bin.
- Ensures input validation for proper binning and error calculations.
- Optimized to handle large event lists efficiently.
"""
function create_lightcurve(eventlist::EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T
    if err_method ∉ (:poisson,)
        throw(ArgumentError("Invalid error method: $err_method. Supported methods: :poisson"))
    end

    isempty(eventlist.times) && throw(ArgumentError("Event list is empty. Cannot create a light curve."))
    bin_size <= 0 && throw(ArgumentError("Bin size must be positive."))
    
    min_time, max_time = extrema(eventlist.times)
    bins = min_time:bin_size:max_time
    
    if length(bins) < 2
        push!(bins, max_time + bin_size)
    end
    
    hist = fit(Histogram, eventlist.times, bins)
    counts = hist.weights 
    
    # to ensure counts vector has correct length
    if length(counts) < length(bins) - 1
        append!(counts, zeros(length(bins) - 1 - length(counts)))
    end
    
    errors = sqrt.(counts) # to aalculate errors based on Poisson statistics

    return LightCurve{T}(bins, counts, errors, err_method)
end