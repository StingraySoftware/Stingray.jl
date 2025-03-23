"""
    LightCurve{T}

A structure representing a light curve, which is a time series of event counts.

## Fields

- `timebins::Vector{T}`: Time bins for the light curve.
- `counts::Vector{T}`: Number of events in each time bin.
- `count_error::Vector{T}`: Error estimate for each bin.
- `err_method::Symbol`: Method used for error estimation (`:poisson`).
"""
struct LightCurve{T}
    timebins::Vector{T}
    counts::Vector{T}
    count_error::Vector{T}
    err_method::Symbol
end

"""
    create_lightcurve(eventlist::EventList{T}, bin_size::T; err_method::Symbol=:poisson)

Generate a light curve from an EventList by binning event times.

## Arguments

- `eventlist::EventList{T}`: The input event list.
- `bin_size::T`: The size of each time bin.
- `err_method::Symbol`: The method for error estimation (`:poisson`).

## Returns

A LightCurve instance containing binned event counts.

## Notes

- The function bins event times into intervals of `bin_size`.
- Errors are estimated using Poisson statistics.
"""
function create_lightcurve(eventlist::EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T
    # Validate err_method
    if err_method != :poisson
        throw(ArgumentError("Invalid error method. Only ':poisson' is supported."))
    end
    min_time = minimum(eventlist.times)
    max_time = maximum(eventlist.times)
    bins = min_time:bin_size:max_time
    hist = fit(Histogram, eventlist.times, bins)
    counts = hist.weights
    errors = sqrt.(counts)
    return LightCurve{T}(bins[1:end-1], counts, errors, err_method)
end
