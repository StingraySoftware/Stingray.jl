module LightCurveModule

import ..Events  # Import EventList from events.jl
export LightCurve, create_lightcurve

"""
    LightCurve{T}

A structure representing a light curve, which is a time series of event counts.

## Fields

- `timebins::Vector{T}`: Time bins for the light curve.
- `counts::Vector{T}`: Number of events in each time bin.
- `count_error::Vector{T}`: Error estimate for each bin.
- `err_method::Symbol`: Method used for error estimation (`:poisson` or `:sqrtN`).
"""
struct LightCurve{T}
    timebins::Vector{T}
    counts::Vector{T}
    count_error::Vector{T}
    err_method::Symbol
end

"""
    create_lightcurve(eventlist::Events.EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T

Generate a light curve from an [`EventList`](@ref) by binning event times.

## Arguments

- `eventlist::Events.EventList{T}`: The input event list.
- `bin_size::T`: The size of each time bin.
- `err_method::Symbol`: The method for error estimation (`:poisson` or `:sqrtN`).

## Returns

A [`LightCurve`](@ref) instance containing binned event counts.

## Notes

- The function bins event times into intervals of `bin_size`.
- Errors are estimated using either Poisson statistics or the square root of counts.
"""
function create_lightcurve(eventlist::Events.EventList{T}, bin_size::T; err_method::Symbol=:poisson) where T
    min_time = minimum(eventlist.times)
    max_time = maximum(eventlist.times)
    
    bins = min_time:bin_size:max_time
    n_bins = length(bins) - 1
    
    counts = zeros(T, n_bins)
    errors = zeros(T, n_bins)
    
    for (i, bin_start) in enumerate(bins[1:end-1])
        bin_end = bins[i + 1]
        events_in_bin = filter(x -> bin_start <= x < bin_end, eventlist.times)
        counts[i] = length(events_in_bin)

        if err_method == :poisson || err_method == :sqrtN
            errors[i] = sqrt(counts[i])
        else
            errors[i] = 0.0
        end
    end
    
    return LightCurve{T}(bins[1:end-1], counts, errors, err_method)
end

end  # module LightCurveModule