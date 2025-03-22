using RecipesBase

module LightCurveModule

import ..Events  # Import EventList from events.jl
using Plots  # Importing the Plots library for graph plotting

export LightCurve, create_lightcurve, plot_lightcurve

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

    # Warning if bin size is too large or too small
    if bin_size > (max_time - min_time) / 10
        @warn "The bin size is relatively large compared to the time range. This may result in a low number of bins."
    elseif bin_size < 0.001
        @warn "The bin size is extremely small. This might lead to excessive error or unnecessary fine-grained bins."
    end

    # Check for valid error method
    if err_method != :poisson && err_method != :sqrtN
        @warn "Unrecognized error method. Using default ':poisson' method."
        err_method = :poisson  # Default to :poisson if invalid
    end

    # Binning events and calculating counts and errors
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

# Plot recipe using RecipesBase
@recipe function f(lightcurve::LightCurve{T}) where T
    # Extract data from the LightCurve object
    timebins = lightcurve.timebins
    counts = lightcurve.counts
    errors = lightcurve.count_error

    # Plot data with error ribbon
    return _plot(
        timebins, counts,
        ribbon=errors,
        label="Event Counts",
        xlabel="Time",
        ylabel="Counts",
        legend=:topright
    )
end

end  # module LightCurveModule
