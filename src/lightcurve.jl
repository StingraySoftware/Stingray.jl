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
    min_time, max_time = extrema(eventlist.times)
    bins = min_time:bin_size:max_time
    hist = fit(Histogram, eventlist.times, bins)
    LightCurve{T}(
        bins[1:end-1] .+ bin_size/2, #timebin
        hist.weights, #count
        sqrt.(hist.weights),#errors
        err_method
    )
end