"""
    LightCurve{T}

A structure representing a light curve, which is a time series of event counts.

## Fields

- `timebins`: Time bins for the light curve.
- `counts`: Number of events in each time bin.
- `count_error`: Error estimate for each bin.
- `err_method`: Method used for error estimation (`:poisson`).
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
    create_lightcurve(eventlist, bin_size; err_method=:poisson)

Generate a light curve by binning event times from an `EventList`.

# Arguments
- `eventlist`: The input event list containing time-stamped events.
- `bin_size`: The size of each time bin.
- `err_method`: The method used for error estimation (currently only `:poisson` is supported).

# Returns
A `LightCurve` object containing the binned event counts and their estimated errors.

# Notes
- Bins event times into intervals of `bin_size` using `fit(Histogram, ...)`.
- Errors are estimated using Poisson statistics (`sqrt(counts)`).
- If no events fall into a bin, the count is zero.
"""
function create_lightcurve(eventlist::EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T
    if err_method ∉ (:poisson,)
        throw(ArgumentError("Invalid error method: $err_method. Supported methods: :poisson"))
    end

    bin_size <= 0 && throw(ArgumentError("Bin size must be positive."))

    min_time, max_time = extrema(eventlist.times)
    bins = min_time:bin_size:max_time

    hist = fit(Histogram, eventlist.times, bins)
    counts = hist.weights

    errors = sqrt.(counts)

    return LightCurve{T}(bins, counts, errors, err_method)
end