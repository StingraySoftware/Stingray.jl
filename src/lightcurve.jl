
"""
Abstract type for all light curve implementations.
"""
abstract type AbstractLightCurve{T} end

"""
    EventProperty{T}

A structure to hold additional event properties beyond time and energy.
"""
struct EventProperty{T}
    name::Symbol
    values::Vector{T}
    unit::String
end

"""
    LightCurveMetadata

A structure containing metadata for light curves.
"""
struct LightCurveMetadata
    telescope::String
    instrument::String
    object::String
    mjdref::Float64
    time_range::Tuple{Float64,Float64}
    bin_size::Float64
    headers::Vector{Dict{String,Any}}
    extra::Dict{String,Any}
end

"""
    LightCurve{T} <: AbstractLightCurve{T}

A structure representing a binned time series with additional properties.
"""
struct LightCurve{T} <: AbstractLightCurve{T}
    timebins::Vector{T}
    bin_edges::Vector{T}
    counts::Vector{Int}
    count_error::Vector{T}
    exposure::Vector{T}
    properties::Vector{EventProperty}
    metadata::LightCurveMetadata
    err_method::Symbol
end

"""
    calculate_errors(counts::Vector{Int}, method::Symbol, exposure::Vector{T}) where T

Calculate statistical uncertainties for count data.
"""
function calculate_errors(
    counts::Vector{Int},
    method::Symbol,
    exposure::Vector{T},
) where {T}
    if method === :poisson
        return convert.(T, sqrt.(counts))
    elseif method === :gaussian
        return convert.(T, sqrt.(counts .+ 1))
    else
        throw(ArgumentError("Unsupported error method: $method. Use :poisson or :gaussian"))
    end
end

"""
    create_lightcurve(
        eventlist::EventList{T}, 
        binsize::Real;
        err_method::Symbol=:poisson,
        tstart::Union{Nothing,Real}=nothing,
        tstop::Union{Nothing,Real}=nothing,
        filters::Dict{Symbol,Any}=Dict{Symbol,Any}()
    ) where T

Create a light curve from an event list with filtering capabilities.
"""
function create_lightcurve(
    eventlist::EventList{T},
    binsize::Real;
    err_method::Symbol = :poisson,
    tstart::Union{Nothing,Real} = nothing,
    tstop::Union{Nothing,Real} = nothing,
    filters::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {T}

    if isempty(eventlist.times)
        throw(ArgumentError("Event list is empty"))
    end

    if binsize <= 0
        throw(ArgumentError("Bin size must be positive"))
    end

    # Initial filtering step
    times = copy(eventlist.times)
    energies = copy(eventlist.energies)

    # Apply time range filter
    start_time = isnothing(tstart) ? minimum(times) : convert(T, tstart)
    stop_time = isnothing(tstop) ? maximum(times) : convert(T, tstop)

    # Filter indices based on all criteria
    valid_indices = findall(t -> start_time ≤ t ≤ stop_time, times)

    # Apply additional filters
    for (key, value) in filters
        if key == :energy
            if value isa Tuple
                energy_indices = findall(e -> value[1] ≤ e < value[2], energies)
                valid_indices = intersect(valid_indices, energy_indices)
            end
        end
    end

    total_events = length(times)
    filtered_events = length(valid_indices)

    #[this below function needs to be discussed properly!]
    # Create bins regardless of whether we have events[i have enter this what if we got unexpectedly allmevents filter out ]
    binsize_t = convert(T, binsize)

    # Make sure we have at least one bin even if start_time equals stop_time
    if start_time == stop_time
        stop_time = start_time + binsize_t
    end

    # Ensure the edges encompass the entire range
    start_bin = floor(start_time / binsize_t) * binsize_t
    num_bins = ceil(Int, (stop_time - start_bin) / binsize_t)
    edges = [start_bin + i * binsize_t for i = 0:num_bins]
    centers = edges[1:end-1] .+ binsize_t / 2

    # Count events in bins
    counts = zeros(Int, length(centers))

    # Only process events if we have any after filtering
    if !isempty(valid_indices)
        filtered_times = times[valid_indices]

        for t in filtered_times
            bin_idx = floor(Int, (t - start_bin) / binsize_t) + 1
            if 1 ≤ bin_idx ≤ length(counts)
                counts[bin_idx] += 1
            end
        end
    end

    # Calculate exposures and errors
    exposure = fill(binsize_t, length(centers))
    errors = calculate_errors(counts, err_method, exposure)

    # Create additional properties
    properties = Vector{EventProperty}()

    # Calculate mean energy per bin if available
    if !isempty(valid_indices) && !isempty(energies)
        filtered_times = times[valid_indices]
        filtered_energies = energies[valid_indices]

        energy_bins = zeros(T, length(centers))
        energy_counts = zeros(Int, length(centers))

        for (t, e) in zip(filtered_times, filtered_energies)
            bin_idx = floor(Int, (t - start_bin) / binsize_t) + 1
            if 1 ≤ bin_idx ≤ length(counts)
                energy_bins[bin_idx] += e
                energy_counts[bin_idx] += 1
            end
        end

        mean_energy = zeros(T, length(centers))
        for i in eachindex(mean_energy)
            mean_energy[i] =
                energy_counts[i] > 0 ? energy_bins[i] / energy_counts[i] : zero(T)
        end

        push!(properties, EventProperty{T}(:mean_energy, mean_energy, "keV"))
    end

    # Create extra metadata with warning if no events remain after filtering
    extra = Dict{String,Any}(
        "filtered_nevents" => filtered_events,
        "total_nevents" => total_events,
        "applied_filters" => filters,
    )

    if filtered_events == 0
        extra["warning"] = "No events remain after filtering"
    end

    # Create metadata
    metadata = LightCurveMetadata(
        get(eventlist.metadata.headers[1], "TELESCOP", ""),
        get(eventlist.metadata.headers[1], "INSTRUME", ""),
        get(eventlist.metadata.headers[1], "OBJECT", ""),
        get(eventlist.metadata.headers[1], "MJDREF", 0.0),
        (start_time, stop_time),
        binsize_t,
        eventlist.metadata.headers,
        extra,
    )

    return LightCurve{T}(
        centers,
        collect(edges),
        counts,
        errors,
        exposure,
        properties,
        metadata,
        err_method,
    )
end

"""
    rebin(lc::LightCurve{T}, new_binsize::Real) where T

Rebin a light curve to a new time resolution.
"""
function rebin(lc::LightCurve{T}, new_binsize::Real) where {T}
    if new_binsize <= lc.metadata.bin_size
        throw(ArgumentError("New bin size must be larger than current bin size"))
    end

    old_binsize = lc.metadata.bin_size
    new_binsize_t = convert(T, new_binsize)

    # Create new bin edges using the same approach as in create_lightcurve
    start_time = lc.metadata.time_range[1]
    stop_time = lc.metadata.time_range[2]

    # Calculate bin edges using the same algorithm as in create_lightcurve
    start_bin = floor(start_time / new_binsize_t) * new_binsize_t
    num_bins = ceil(Int, (stop_time - start_bin) / new_binsize_t)
    new_edges = [start_bin + i * new_binsize_t for i = 0:num_bins]
    new_centers = new_edges[1:end-1] .+ new_binsize_t / 2

    # Rebin counts
    new_counts = zeros(Int, length(new_centers))

    for (i, time) in enumerate(lc.timebins)
        if lc.counts[i] > 0  # Only process bins with counts
            bin_idx = floor(Int, (time - start_bin) / new_binsize_t) + 1
            if 1 ≤ bin_idx ≤ length(new_counts)
                new_counts[bin_idx] += lc.counts[i]
            end
        end
    end

    # Calculate new exposures and errors
    new_exposure = fill(new_binsize_t, length(new_centers))
    new_errors = calculate_errors(new_counts, lc.err_method, new_exposure)

    # Rebin properties
    new_properties = Vector{EventProperty}()
    for prop in lc.properties
        new_values = zeros(T, length(new_centers))
        counts = zeros(Int, length(new_centers))

        for (i, val) in enumerate(prop.values)
            if lc.counts[i] > 0  # Only process bins with counts
                bin_idx = floor(Int, (lc.timebins[i] - start_bin) / new_binsize_t) + 1
                if 1 ≤ bin_idx ≤ length(new_values)
                    new_values[bin_idx] += val * lc.counts[i]
                    counts[bin_idx] += lc.counts[i]
                end
            end
        end

        # Calculate weighted average
        for i in eachindex(new_values)
            new_values[i] = counts[i] > 0 ? new_values[i] / counts[i] : zero(T)
        end

        push!(new_properties, EventProperty(prop.name, new_values, prop.unit))
    end

    # Update metadata
    new_metadata = LightCurveMetadata(
        lc.metadata.telescope,
        lc.metadata.instrument,
        lc.metadata.object,
        lc.metadata.mjdref,
        lc.metadata.time_range,
        new_binsize_t,
        lc.metadata.headers,
        merge(lc.metadata.extra, Dict{String,Any}("original_binsize" => old_binsize)),
    )

    return LightCurve{T}(
        new_centers,
        collect(new_edges),
        new_counts,
        new_errors,
        new_exposure,
        new_properties,
        new_metadata,
        lc.err_method,
    )
end

# Basic array interface methods
Base.length(lc::LightCurve) = length(lc.counts)
Base.size(lc::LightCurve) = (length(lc),)
Base.getindex(lc::LightCurve, i) = (lc.timebins[i], lc.counts[i])
