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
    calculate_errors(counts::Vector{Int}, method::Symbol, exposure::Vector{T}; 
                    gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T

Calculate statistical uncertainties for count data using vectorized operations.
"""
function calculate_errors(counts::Vector{Int}, method::Symbol, exposure::Vector{T}; 
                         gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T
    if method === :poisson
        # Vectorized Poisson errors: σ = sqrt(N), use sqrt(N + 1) when N = 0
        return convert.(T, @. sqrt(max(counts, 1)))
    elseif method === :gaussian
        if isnothing(gaussian_errors)
            throw(ArgumentError("Gaussian errors must be provided by user when using :gaussian method"))
        end
        if length(gaussian_errors) != length(counts)
            throw(ArgumentError("Length of gaussian_errors must match length of counts"))
        end
        return gaussian_errors
    else
        throw(ArgumentError("Unsupported error method: $method. Use :poisson or :gaussian"))
    end
end

"""
    validate_lightcurve_inputs(eventlist, binsize, err_method, gaussian_errors)

Validate all inputs for light curve creation before processing.
"""
function validate_lightcurve_inputs(eventlist, binsize, err_method, gaussian_errors)
    # Check event list
    if isempty(eventlist.times)
        throw(ArgumentError("Event list is empty"))
    end
    
    # Check bin size
    if binsize <= 0
        throw(ArgumentError("Bin size must be positive"))
    end
    
    # Check error method
    if !(err_method in [:poisson, :gaussian])
        throw(ArgumentError("Unsupported error method: $err_method. Use :poisson or :gaussian"))
    end
    
    # Check Gaussian errors if needed
    if err_method === :gaussian
        if isnothing(gaussian_errors)
            throw(ArgumentError("Gaussian errors must be provided when using :gaussian method"))
        end
        # Note: Length validation will happen after filtering, not here
    end
end

"""
    apply_event_filters(times::Vector{T}, energies::Union{Nothing,Vector{T}}, 
                       tstart::Union{Nothing,Real}, tstop::Union{Nothing,Real},
                       energy_filter::Union{Nothing,Tuple{Real,Real}}) where T

Apply time and energy filters to event data.
Returns filtered times and energies.
"""
function apply_event_filters(times::Vector{T}, energies::Union{Nothing,Vector{T}}, 
                            tstart::Union{Nothing,Real}, tstop::Union{Nothing,Real},
                            energy_filter::Union{Nothing,Tuple{Real,Real}}) where T
    
    filtered_times = times
    filtered_energies = energies
    
    # Apply energy filter first if specified
    if !isnothing(energy_filter) && !isnothing(energies)
        emin, emax = energy_filter
        energy_mask = @. (energies >= emin) & (energies < emax)
        filtered_times = times[energy_mask]
        filtered_energies = energies[energy_mask]
        
        if isempty(filtered_times)
            throw(ArgumentError("No events remain after energy filtering"))
        end
        @info "Applied energy filter [$emin, $emax) keV: $(length(filtered_times)) events remain"
    end
    
    # Determine time range
    start_time = isnothing(tstart) ? minimum(filtered_times) : convert(T, tstart)
    stop_time = isnothing(tstop) ? maximum(filtered_times) : convert(T, tstop)
    
    # Apply time filter if needed
    if start_time != minimum(filtered_times) || stop_time != maximum(filtered_times)
        time_mask = @. (filtered_times >= start_time) & (filtered_times <= stop_time)
        filtered_times = filtered_times[time_mask]
        if !isnothing(filtered_energies)
            filtered_energies = filtered_energies[time_mask]
        end
        
        if isempty(filtered_times)
            throw(ArgumentError("No events remain after time filtering"))
        end
        @info "Applied time filter [$start_time, $stop_time]: $(length(filtered_times)) events remain"
    end
    
    return filtered_times, filtered_energies, start_time, stop_time
end

"""
    create_time_bins(start_time::T, stop_time::T, binsize::T) where T

Create time bin edges and centers for the light curve.
"""
function create_time_bins(start_time::T, stop_time::T, binsize::T) where T
    # Ensure we cover the full range including the endpoint
    start_bin = floor(start_time / binsize) * binsize
    
    # Calculate number of bins to ensure we cover stop_time
    time_span = stop_time - start_bin
    num_bins = max(1, ceil(Int, time_span / binsize))
    
    # Adjust if the calculated end would be less than stop_time
    while start_bin + num_bins * binsize < stop_time
        num_bins += 1
    end
    
    # Create bin edges and centers efficiently
    edges = [start_bin + i * binsize for i in 0:num_bins]
    centers = [start_bin + (i + 0.5) * binsize for i in 0:(num_bins-1)]
    
    return edges, centers
end

"""
    bin_events(times::Vector{T}, bin_edges::Vector{T}) where T

Bin event times into histogram counts.
"""
function bin_events(times::Vector{T}, bin_edges::Vector{T}) where T
    # Use StatsBase for fast, memory-efficient binning
    hist = fit(Histogram, times, bin_edges)
    return Vector{Int}(hist.weights)
end

"""
    calculate_additional_properties(times::Vector{T}, energies::Union{Nothing,Vector{U}}, 
                                   bin_edges::Vector{T}, bin_centers::Vector{T}) where {T,U}

Calculate additional properties like mean energy per bin.
handles type mismatches between time and energy vectors.
"""
function calculate_additional_properties(times::Vector{T}, energies::Union{Nothing,Vector{U}}, 
                                        bin_edges::Vector{T}, bin_centers::Vector{T}) where {T,U}
    properties = Vector{EventProperty}()
    
    # Calculate mean energy per bin if available
    if !isnothing(energies) && !isempty(energies) && length(bin_centers) > 0
        start_bin = bin_edges[1]
        
        # Handle case where there's only one bin center
        if length(bin_centers) == 1
            binsize = length(bin_edges) > 1 ? bin_edges[2] - bin_edges[1] : T(1)
        else
            binsize = bin_centers[2] - bin_centers[1]  # Assuming uniform bins
        end
        
        # Use efficient binning for energies
        energy_sums = zeros(T, length(bin_centers))
        energy_counts = zeros(Int, length(bin_centers))
        
        # Vectorized binning for energies
        for (t, e) in zip(times, energies)
            bin_idx = floor(Int, (t - start_bin) / binsize) + 1
            if 1 ≤ bin_idx ≤ length(bin_centers)
                energy_sums[bin_idx] += T(e)  # Convert energy to time type
                energy_counts[bin_idx] += 1
            end
        end
        
        # Calculate mean energies using vectorized operations
        mean_energy = @. ifelse(energy_counts > 0, energy_sums / energy_counts, zero(T))
        push!(properties, EventProperty{T}(:mean_energy, mean_energy, "keV"))
    end
    
    return properties
end

"""
    extract_metadata(eventlist, start_time, stop_time, binsize, filtered_times, energy_filter)

Extract and create metadata for the light curve.
"""
function extract_metadata(eventlist, start_time, stop_time, binsize, filtered_times, energy_filter)
    first_header = isempty(eventlist.metadata.headers) ? Dict{String,Any}() : eventlist.metadata.headers[1]
    
    return LightCurveMetadata(
        get(first_header, "TELESCOP", ""),
        get(first_header, "INSTRUME", ""),
        get(first_header, "OBJECT", ""),
        get(first_header, "MJDREF", 0.0),
        (Float64(start_time), Float64(stop_time)),
        Float64(binsize),
        eventlist.metadata.headers,
        Dict{String,Any}(
            "filtered_nevents" => length(filtered_times),
            "total_nevents" => length(eventlist.times),
            "energy_filter" => energy_filter
        )
    )
end

"""
    create_lightcurve(
        eventlist::EventList{T}, 
        binsize::Real;
        err_method::Symbol=:poisson,
        gaussian_errors::Union{Nothing,Vector{T}}=nothing,
        tstart::Union{Nothing,Real}=nothing,
        tstop::Union{Nothing,Real}=nothing,
        energy_filter::Union{Nothing,Tuple{Real,Real}}=nothing,
        event_filter::Union{Nothing,Function}=nothing
    ) where T

Create a light curve from an event list with enhanced performance and filtering.

# Arguments
- `eventlist`: The input event list
- `binsize`: Time bin size
- `err_method`: Error calculation method (:poisson or :gaussian)
- `gaussian_errors`: User-provided Gaussian errors (required if err_method=:gaussian)
- `tstart`, `tstop`: Time range limits
- `energy_filter`: Energy range as (emin, emax) tuple
- `event_filter`: Optional function to filter events, should return boolean mask
"""
function create_lightcurve(
    eventlist::EventList{T}, 
    binsize::Real;
    err_method::Symbol=:poisson,
    gaussian_errors::Union{Nothing,Vector{T}}=nothing,
    tstart::Union{Nothing,Real}=nothing,
    tstop::Union{Nothing,Real}=nothing,
    energy_filter::Union{Nothing,Tuple{Real,Real}}=nothing,
    event_filter::Union{Nothing,Function}=nothing
) where T
    
    # Validate all inputs first
    validate_lightcurve_inputs(eventlist, binsize, err_method, gaussian_errors)
    
    binsize_t = convert(T, binsize)
    
    # Get initial data references
    times = eventlist.times
    energies = eventlist.energies
    
    # Apply custom event filter if provided
    if !isnothing(event_filter)
        filter_mask = event_filter(eventlist)
        if !isa(filter_mask, AbstractVector{Bool})
            throw(ArgumentError("Event filter function must return a boolean vector"))
        end
        if length(filter_mask) != length(times)
            throw(ArgumentError("Event filter mask length must match number of events"))
        end
        
        times = times[filter_mask]
        if !isnothing(energies)
            energies = energies[filter_mask]
        end
        
        if isempty(times)
            throw(ArgumentError("No events remain after custom filtering"))
        end
        @info "Applied custom filter: $(length(times)) events remain"
    end
    
    # Apply standard filters
    filtered_times, filtered_energies, start_time, stop_time = apply_event_filters(
        times, energies, tstart, tstop, energy_filter
    )
    
    # Create time bins
    bin_edges, bin_centers = create_time_bins(start_time, stop_time, binsize_t)
    
    # Bin the events
    counts = bin_events(filtered_times, bin_edges)
    
    @info "Created light curve: $(length(bin_centers)) bins, bin size = $(binsize_t) s"
    
    # Now validate gaussian_errors length if needed
    if err_method === :gaussian && !isnothing(gaussian_errors)
        if length(gaussian_errors) != length(counts)
            throw(ArgumentError("Length of gaussian_errors ($(length(gaussian_errors))) must match number of bins ($(length(counts)))"))
        end
    end
    
    # Calculate exposures and errors
    exposure = fill(binsize_t, length(bin_centers))
    errors = calculate_errors(counts, err_method, exposure; gaussian_errors=gaussian_errors)
    
    # Calculate additional properties
    properties = calculate_additional_properties(filtered_times, filtered_energies, bin_edges, bin_centers)
    
    # Extract metadata
    metadata = extract_metadata(eventlist, start_time, stop_time, binsize_t, filtered_times, energy_filter)
    
    return LightCurve{T}(
        bin_centers,
        bin_edges,
        counts,
        errors,
        exposure,
        properties,
        metadata,
        err_method
    )
end

"""
    rebin(lc::LightCurve{T}, new_binsize::Real; 
          gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T

Rebin a light curve to a new time resolution with enhanced performance.
"""
function rebin(lc::LightCurve{T}, new_binsize::Real; 
               gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T
    if new_binsize <= lc.metadata.bin_size
        throw(ArgumentError("New bin size must be larger than current bin size"))
    end
    
    old_binsize = T(lc.metadata.bin_size)
    new_binsize_t = convert(T, new_binsize)
    
    # Create new bin edges using the same approach as in create_lightcurve
    start_time = T(lc.metadata.time_range[1])
    stop_time = T(lc.metadata.time_range[2])
    
    # Calculate bin edges using efficient algorithm
    start_bin = floor(start_time / new_binsize_t) * new_binsize_t
    time_span = stop_time - start_bin
    num_bins = max(1, ceil(Int, time_span / new_binsize_t))
    
    # Ensure we cover the full range
    while start_bin + num_bins * new_binsize_t < stop_time
        num_bins += 1
    end
    
    new_edges = [start_bin + i * new_binsize_t for i in 0:num_bins]
    new_centers = [start_bin + (i + 0.5) * new_binsize_t for i in 0:(num_bins-1)]
    
    # Rebin counts using vectorized operations where possible
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
    
    # Handle error propagation based on original method
    if lc.err_method === :gaussian && isnothing(gaussian_errors)
        throw(ArgumentError("Gaussian errors must be provided when rebinning a light curve with Gaussian errors"))
    end
    
    new_errors = calculate_errors(new_counts, lc.err_method, new_exposure; gaussian_errors=gaussian_errors)
    
    # Rebin properties using weighted averaging
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
        
        # Calculate weighted average using vectorized operations
        new_values = @. ifelse(counts > 0, new_values / counts, zero(T))
        
        push!(new_properties, EventProperty(prop.name, new_values, prop.unit))
    end
    
    # Update metadata
    new_metadata = LightCurveMetadata(
        lc.metadata.telescope,
        lc.metadata.instrument,
        lc.metadata.object,
        lc.metadata.mjdref,
        lc.metadata.time_range,
        Float64(new_binsize_t),
        lc.metadata.headers,
        merge(
            lc.metadata.extra,
            Dict{String,Any}("original_binsize" => Float64(old_binsize))
        )
    )
    
    return LightCurve{T}(
        new_centers,
        new_edges,
        new_counts,
        new_errors,
        new_exposure,
        new_properties,
        new_metadata,
        lc.err_method
    )
end

# Basic array interface methods
Base.length(lc::LightCurve) = length(lc.counts)
Base.size(lc::LightCurve) = (length(lc),)
Base.getindex(lc::LightCurve, i) = (lc.timebins[i], lc.counts[i])