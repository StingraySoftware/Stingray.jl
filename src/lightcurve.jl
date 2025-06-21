"""
    abstract type AbstractLightCurve{T} end

# Examples
```julia
# All light curve types inherit from this
struct MyLightCurve{T} <: AbstractLightCurve{T}
    # ... fields
end
```
"""
abstract type AbstractLightCurve{T} end

"""
    EventProperty{T}

A structure to hold additional event properties beyond time and energy.

This structure stores computed properties that are calculated per time bin,
such as mean energy, hardness ratios, or other derived quantities.

# Fields
- `name::Symbol`: Name identifier for the property
- `values::Vector{T}`: Property values for each time bin
- `unit::String`: Physical units of the property values

# Examples
```julia
# Create a property for mean energy per bin
mean_energy = EventProperty{Float64}(:mean_energy, [1.2, 1.5, 1.8], "keV")

# Create a property for count rates
count_rate = EventProperty{Float64}(:rate, [10.5, 12.1, 9.8], "counts/s")
```
"""
struct EventProperty{T}
    "Name identifier for the property"
    name::Symbol
    "Property values for each time bin"
    values::Vector{T}
    "Physical units of the property values"
    unit::String
end

"""
    LightCurveMetadata

A structure containing comprehensive metadata for light curves.

This structure stores all relevant metadata about the light curve creation,
including source information, timing parameters, and processing history.

# Fields
- `telescope::String`: Name of the telescope/mission
- `instrument::String`: Name of the instrument
- `object::String`: Name of the observed object
- `mjdref::Float64`: Modified Julian Date reference time
- `time_range::Tuple{Float64,Float64}`: Start and stop times of the light curve
- `bin_size::Float64`: Time bin size in seconds
- `headers::Vector{Dict{String,Any}}`: Original FITS headers
- `extra::Dict{String,Any}`: Additional metadata and processing information

# Examples
```julia
# Metadata is typically created automatically
lc = create_lightcurve(eventlist, 1.0)
println(lc.metadata.telescope)    # "NICER"
println(lc.metadata.bin_size)     # 1.0
println(lc.metadata.time_range)   # (1000.0, 2000.0)
```
"""
struct LightCurveMetadata
    "Name of the telescope/mission"
    telescope::String
    "Name of the instrument"
    instrument::String
    "Name of the observed object"
    object::String
    "Modified Julian Date reference time"
    mjdref::Float64
    "Start and stop times of the light curve"
    time_range::Tuple{Float64,Float64}
    "Time bin size in seconds"
    bin_size::Float64
    "Original FITS headers"
    headers::Vector{Dict{String,Any}}
    "Additional metadata and processing information"
    extra::Dict{String,Any}
end

"""
    LightCurve{T} <: AbstractLightCurve{T}

A structure representing a binned time series with additional properties.

This is the main light curve structure that holds binned photon count data
along with statistical uncertainties, exposure times, and derived properties.

# Fields
- `timebins::Vector{T}`: Time bin centers
- `bin_edges::Vector{T}`: Time bin edges (length = timebins + 1)
- `counts::Vector{Int}`: Photon counts in each bin
- `count_error::Vector{T}`: Statistical uncertainties on counts
- `exposure::Vector{T}`: Exposure time for each bin
- `properties::Vector{EventProperty}`: Additional computed properties
- `metadata::LightCurveMetadata`: Comprehensive metadata
- `err_method::Symbol`: Error calculation method used (:poisson or :gaussian)

# Examples
```julia
# Create from event list
ev = readevents("events.fits")
lc = create_lightcurve(ev, 1.0)  # 1-second bins

# Access data
println("Time bins: ", lc.timebins[1:5])
println("Counts: ", lc.counts[1:5])
println("Errors: ", lc.count_error[1:5])

# Basic operations
println("Total counts: ", sum(lc.counts))
println("Mean count rate: ", mean(lc.counts ./ lc.exposure))
```

# Interface
- `length(lc)`: Number of time bins
- `lc[i]`: Get (time, counts) tuple for bin i
- Supports array-like indexing and iteration

See also [`create_lightcurve`](@ref), [`rebin`](@ref).
"""
struct LightCurve{T} <: AbstractLightCurve{T}
    "Time bin centers"
    timebins::Vector{T}
    "Time bin edges (length = timebins + 1)"
    bin_edges::Vector{T}
    "Photon counts in each bin"
    counts::Vector{Int}
    "Statistical uncertainties on counts"
    count_error::Vector{T}
    "Exposure time for each bin"
    exposure::Vector{T}
    "Additional computed properties"
    properties::Vector{EventProperty{T}}
    "Comprehensive metadata"
    metadata::LightCurveMetadata
    "Error calculation method used (:poisson or :gaussian)"
    err_method::Symbol
end

"""
    calculate_errors(counts::Vector{Int}, method::Symbol, exposure::Vector{T}; 
                    gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T

Calculate statistical uncertainties for count data using vectorized operations.

This function computes appropriate statistical uncertainties based on the
specified method, with optimized vectorized implementations for performance.

# Arguments
- `counts::Vector{Int}`: Photon counts in each bin
- `method::Symbol`: Error calculation method (:poisson or :gaussian)
- `exposure::Vector{T}`: Exposure times (currently unused, kept for interface compatibility)

# Keyword Arguments
- `gaussian_errors::Union{Nothing,Vector{T}}`: User-provided errors for :gaussian method

# Returns
`Vector{T}`: Statistical uncertainties for each bin

# Methods
- `:poisson`: Uses Poisson statistics (σ = √N, with σ = 1 for N = 0)
- `:gaussian`: Uses user-provided Gaussian errors

# Examples
```julia
counts = [10, 25, 5, 0, 15]
exposures = fill(1.0, 5)

# Poisson errors
errors = calculate_errors(counts, :poisson, exposures)
# Result: [3.16, 5.0, 2.24, 1.0, 3.87]

# Gaussian errors
gaussian_errs = [0.5, 0.8, 0.3, 0.1, 0.6]
errors = calculate_errors(counts, :gaussian, exposures; gaussian_errors=gaussian_errs)
# Result: [0.5, 0.8, 0.3, 0.1, 0.6]
```

# Throws
- `ArgumentError`: If method is not :poisson or :gaussian
- `ArgumentError`: If :gaussian method used without providing gaussian_errors
- `ArgumentError`: If gaussian_errors length doesn't match counts length

# Implementation Notes
- Uses vectorized operations with `@.` macro for performance
- Handles zero counts case for Poisson statistics
- Type-stable with explicit type conversions
"""
function calculate_errors(
    counts::Vector{Int},
    method::Symbol,
    exposure::Vector{T};
    gaussian_errors::Union{Nothing,Vector{T}} = nothing,
) where {T}
    if method === :poisson
        # Vectorized Poisson errors: σ = sqrt(N), use sqrt(N + 1) when N = 0
        return convert.(T, @. sqrt(max(counts, 1)))
    elseif method === :gaussian
        if isnothing(gaussian_errors)
            throw(
                ArgumentError(
                    "Gaussian errors must be provided by user when using :gaussian method",
                ),
            )
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

This function performs comprehensive validation of all input parameters
to ensure they are suitable for light curve creation, providing clear
error messages for common issues.

# Arguments
- `eventlist`: EventList structure containing photon arrival times
- `binsize`: Time bin size (must be positive)
- `err_method`: Error calculation method (:poisson or :gaussian)
- `gaussian_errors`: User-provided errors (required for :gaussian method)

# Throws
- `ArgumentError`: If event list is empty
- `ArgumentError`: If bin size is not positive
- `ArgumentError`: If error method is not supported
- `ArgumentError`: If :gaussian method used without providing errors

# Examples
```julia
# This function is called internally by create_lightcurve
# Manual validation for custom workflows:
validate_lightcurve_inputs(ev, 1.0, :poisson, nothing)  # OK
validate_lightcurve_inputs(ev, -1.0, :poisson, nothing) # Throws error
```

# Implementation Notes
- Validates inputs early to provide clear error messages
- Length validation for gaussian_errors happens after filtering
- Separated for testability and modularity
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
        throw(
            ArgumentError(
                "Unsupported error method: $err_method. Use :poisson or :gaussian",
            ),
        )
    end

    # Check Gaussian errors if needed - but don't validate length here
    # Length validation will happen after binning when we know the actual number of bins
    if err_method === :gaussian
        if isnothing(gaussian_errors)
            throw(
                ArgumentError(
                    "Gaussian errors must be provided when using :gaussian method",
                ),
            )
        end
    end
end

"""
    apply_event_filters(times::Vector{T}, energies::Union{Nothing,Vector{T}}, 
                       tstart::Union{Nothing,Real}, tstop::Union{Nothing,Real},
                       energy_filter::Union{Nothing,Tuple{Real,Real}}) where T

Apply time and energy filters to event data with comprehensive validation.

This function applies filtering operations to event data, supporting both
time range and energy range filtering with comprehensive validation and
logging of the filtering process. Filters are applied in the optimal order
to maximize efficiency.

# Arguments
- `times::Vector{T}`: Event arrival times
- `energies::Union{Nothing,Vector{T}}`: Event energies (or nothing)
- `tstart::Union{Nothing,Real}`: Start time for filtering (or nothing for no limit)
- `tstop::Union{Nothing,Real}`: Stop time for filtering (or nothing for no limit)
- `energy_filter::Union{Nothing,Tuple{Real,Real}}`: Energy range as (emin, emax)

# Returns
`Tuple{Vector{T}, Union{Nothing,Vector{T}}, T, T}`: 
- Filtered times
- Filtered energies (or nothing if no energy data)
- Final start time used
- Final stop time used

# Examples
```julia
# Filter by energy only
filtered_times, filtered_energies, start_t, stop_t = apply_event_filters(
    times, energies, nothing, nothing, (0.5, 10.0)
)

# Filter by time only  
filtered_times, filtered_energies, start_t, stop_t = apply_event_filters(
    times, energies, 1000.0, 2000.0, nothing
)

# Filter by both time and energy
filtered_times, filtered_energies, start_t, stop_t = apply_event_filters(
    times, energies, 1000.0, 2000.0, (0.5, 10.0)
)
```

# Throws
- `ArgumentError`: If no events remain after energy filtering
- `ArgumentError`: If no events remain after time filtering

# Implementation Notes
- Applies energy filter first, then time filter for optimal performance
- Uses vectorized boolean operations for efficiency
- Provides informative logging of filtering results
- Handles the case where energies might be nothing
- Automatically determines time range if not specified
- Type-stable implementation with proper type inference
"""
function apply_event_filters(
    times::TimeType,
    energies::Union{Nothing,TimeType},
    tstart::Union{Nothing,Real},
    tstop::Union{Nothing,Real},
    energy_filter::Union{Nothing,Tuple{Real,Real}},
) where {TimeType<:AbstractVector}

    T = eltype(TimeType)

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

Create uniform time bin edges and centers for light curve binning.

This function creates a uniform time grid that covers the specified time range
with the given bin size, ensuring complete coverage of the data range with
proper alignment to bin boundaries.

# Arguments
- `start_time::T`: Start time of the data range
- `stop_time::T`: Stop time of the data range  
- `binsize::T`: Size of each time bin

# Returns
`Tuple{Vector{T}, Vector{T}}`: (bin_edges, bin_centers)
- `bin_edges`: Array of bin boundaries (length = n_bins + 1)
- `bin_centers`: Array of bin center times (length = n_bins)

# Examples
```julia
# Create 1-second bins from 1000 to 1100 seconds
edges, centers = create_time_bins(1000.0, 1100.0, 1.0)
println(length(edges))    # 101 (100 bins + 1)
println(length(centers))  # 100

# First and last bins
println(edges[1])         # 1000.0
println(edges[end])       # 1100.0
println(centers[1])       # 1000.5
println(centers[end])     # 1099.5

# Sub-second binning
edges, centers = create_time_bins(0.0, 10.0, 0.1)
println(length(centers))  # 100 (0.1-second bins)
```

# Implementation Notes
- Aligns start_bin to bin_size boundaries for consistent binning
- Ensures complete coverage of the stop_time
- Uses efficient list comprehensions for bin creation
- Handles edge cases where time span is less than one bin
- Centers are calculated as bin_start + 0.5 * bin_size
- Type-stable implementation preserving input types

# Algorithm Details
1. Aligns starting bin edge to multiple of bin_size before start_time
2. Calculates minimum number of bins needed to cover full range
3. Adds extra bin if needed to ensure stop_time is included
4. Generates bin edges and centers using vectorized operations
"""
function create_time_bins(start_time::T, stop_time::T, binsize::T) where {T}
    # Ensure we cover the full range including the endpoint
    start_bin = floor(start_time / binsize) * binsize

    # Calculate number of bins to ensure we cover stop_time
    # Add a small epsilon to ensure the stop_time is included
    time_span = stop_time - start_bin
    num_bins = max(1, ceil(Int, time_span / binsize))

    # Ensure the last bin includes stop_time by checking if we need an extra bin
    if start_bin + num_bins * binsize <= stop_time
        num_bins += 1
    end

    # Create bin edges and centers efficiently
    edges = [start_bin + i * binsize for i = 0:num_bins]
    centers = [start_bin + (i + 0.5) * binsize for i = 0:(num_bins-1)]

    return edges, centers
end

"""
    bin_events(times::Vector{T}, bin_edges::Vector{T}) where T

Bin event arrival times into histogram counts using optimized algorithms.

This function efficiently bins photon arrival times into a histogram using
optimized algorithms from StatsBase.jl for maximum performance. The last
bin edge is made inclusive to ensure all events are captured.

# Arguments
- `times::Vector{T}`: Event arrival times (need not be sorted)
- `bin_edges::Vector{T}`: Time bin edges (length = n_bins + 1)

# Returns
`Vector{Int}`: Photon counts in each bin

# Examples
```julia
times = [1.1, 1.3, 1.7, 2.2, 2.8, 3.1]
edges = [1.0, 2.0, 3.0, 4.0]  # 3 bins: [1,2), [2,3), [3,4)

counts = bin_events(times, edges)
# Result: [3, 2, 1] (3 events in first bin, 2 in second, 1 in third)

# Edge cases
empty_times = Float64[]
counts = bin_events(empty_times, edges)
# Result: [0, 0, 0] (all bins empty)

# Single event at bin edge
times = [2.0]  # Exactly on bin boundary
counts = bin_events(times, edges)
# Result: [0, 1, 0] (event goes in second bin)
```

# Throws
- `ArgumentError`: If fewer than 2 bin edges provided

# Implementation Notes
- Uses StatsBase.Histogram for optimized binning
- Makes the rightmost bin edge inclusive by adding small epsilon
- Handles edge cases (empty arrays, single events) gracefully
- Results are deterministic and reproducible
- Memory-efficient implementation
- Type-stable with explicit conversion to Vector{Int}

# Performance
- O(n log m) complexity where n = events, m = bins (for unsorted data)
- O(n + m) complexity for pre-sorted data
- Vectorized operations minimize memory allocation
- Efficient for large event lists and many bins

# Binning Rules
- All bins except the last are half-open intervals [a, b)
- The last bin is closed interval [a, b] to capture boundary events
- Events exactly on interior boundaries belong to the right bin
"""
function bin_events(
    times::TimeType,
    bin_edges::Vector{T},
) where {TimeType<:AbstractVector,T}
    if length(bin_edges) < 2
        throw(ArgumentError("Need at least 2 bin edges"))
    end

    # Use StatsBase histogram but ensure the rightmost edge is inclusive
    # by slightly expanding the last edge
    adjusted_edges = copy(bin_edges)
    if length(adjusted_edges) > 1
        # Add small epsilon to make the last bin inclusive of the right edge
        adjusted_edges[end] = nextfloat(adjusted_edges[end])
    end

    hist = fit(Histogram, times, adjusted_edges)
    return Vector{Int}(hist.weights)
end

"""
    calculate_additional_properties(times::Vector{T}, energies::Union{Nothing,Vector{U}}, 
                                   bin_edges::Vector{T}, bin_centers::Vector{T}) where {T,U}

Calculate derived properties per time bin from event data.

This function computes derived properties for each time bin, currently
focusing on energy-related statistics but extensible to other properties.
All calculations use efficient vectorized operations where possible.

# Arguments
- `times::Vector{T}`: Event arrival times
- `energies::Union{Nothing,Vector{U}}`: Event energies (or nothing if unavailable)
- `bin_edges::Vector{T}`: Time bin edges
- `bin_centers::Vector{T}`: Time bin centers

# Returns
`Vector{EventProperty}`: Vector of computed properties

# Properties Calculated
- `mean_energy`: Average photon energy per time bin (if energy data available)

# Examples
```julia
# With energy data
times = [1.1, 1.3, 2.2, 2.8]
energies = [1.5, 2.0, 1.8, 2.2]
edges = [1.0, 2.0, 3.0]
centers = [1.5, 2.5]

props = calculate_additional_properties(times, energies, edges, centers)
# Result: [EventProperty(:mean_energy, [1.75, 2.0], "keV")]

# Without energy data
props = calculate_additional_properties(times, nothing, edges, centers)
# Result: [] (empty vector)

# Empty bins handled gracefully
times = [1.1]
energies = [1.5]
edges = [1.0, 2.0, 3.0, 4.0]  # 3 bins, only first has events
centers = [1.5, 2.5, 3.5]

props = calculate_additional_properties(times, energies, edges, centers)
# Result: [EventProperty(:mean_energy, [1.5, 0.0, 0.0], "keV")]
```

# Implementation Notes
- Handles type mismatches between time and energy vectors gracefully
- Uses efficient vectorized operations and pre-allocated arrays
- Gracefully handles edge cases (empty bins, missing data, single bins)
- Extensible design for adding new properties
- Type-stable with explicit type conversions
- Zero values assigned to bins with no events

# Performance
- O(n) complexity where n = number of events
- Minimal memory allocation through pre-allocated arrays
- Vectorized operations for arithmetic computations
- Efficient indexing for bin assignment

# Future Extensions
This function can be extended to calculate additional properties such as:
- Hardness ratios between energy bands
- Energy spread (standard deviation) per bin
- Spectral indices or colors
- Custom user-defined properties via callback functions

# Binning Algorithm
Events are assigned to bins using:
```julia
bin_idx = floor(Int, (time - start_bin) / binsize) + 1
```
This ensures consistent binning with `create_time_bins` and `bin_events`.
"""
function calculate_additional_properties(
    times::TimeType,
    energies::Union{Nothing,EnergyType},
    bin_edges::Vector{T},
    bin_centers::Vector{T},
) where {TimeType<:AbstractVector,EnergyType<:AbstractVector,T}
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
            if 1 <= bin_idx ≤ length(bin_centers)
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

Extract and organize metadata from event list and processing parameters.

This function creates comprehensive metadata for the light curve by extracting
information from the original event list headers and combining it with
processing parameters and filtering information.

# Arguments
- `eventlist`: Input EventList structure
- `start_time`: Final start time after filtering
- `stop_time`: Final stop time after filtering  
- `binsize`: Time bin size used
- `filtered_times`: Final filtered event times
- `energy_filter`: Energy filter applied (or nothing)

# Returns
`LightCurveMetadata`: Complete metadata structure

# Implementation Notes
- Preserves ALL original FITS headers for full traceability
- Handles various header formats (FITSHeader, Vector, Dict)
- Extracts common astronomical fields with sensible defaults
- Records complete processing history in extra metadata
- Maintains backward compatibility with different input formats

# Metadata Fields Extracted
- Standard FITS keywords: TELESCOP, INSTRUME, OBJECT, MJDREF
- Processing information: event counts, filters applied
- Timing information: binsize, time range
- Complete header preservation for reference

# Examples
```julia
# Typically called internally by create_lightcurve
metadata = extract_metadata(ev, 1000.0, 2000.0, 1.0, filtered_times, (0.5, 10.0))
println(metadata.telescope)  # "NICER"
println(metadata.extra["filtered_nevents"])  # 15000
```
"""

function extract_metadata(
    eventlist,
    start_time,
    stop_time,
    binsize,
    filtered_times,
    energy_filter,
)
    # Convert headers to the expected format - preserve ALL original metadata
    headers = if eventlist.meta.headers isa FITSIO.FITSHeader
        [Dict{String,Any}(pairs(eventlist.meta.headers))]
    elseif eventlist.meta.headers isa Vector
        eventlist.meta.headers
    elseif eventlist.meta.headers isa Dict
        [eventlist.meta.headers]
    else
        [Dict{String,Any}()]
    end

    first_header = isempty(headers) ? Dict{String,Any}() : headers[1]

    # Extract common astronomical fields with defaults, but don't force specific values
    telescope = get(first_header, "TELESCOP", get(first_header, "TELESCOPE", ""))
    instrument = get(first_header, "INSTRUME", get(first_header, "INSTRUMENT", ""))
    object = get(first_header, "OBJECT", get(first_header, "TARGET", ""))
    mjdref = get(first_header, "MJDREF", get(first_header, "MJDREFI", 0.0))

    # Create comprehensive extra metadata including processing info
    extra_metadata = Dict{String,Any}(
        "filtered_nevents" => length(filtered_times),
        "total_nevents" => length(eventlist.times),
        "energy_filter" => energy_filter,
        "binning_method" => "histogram",
    )

    # Add any additional metadata from the eventlist that's not in headers
    if hasfield(typeof(eventlist.meta), :extra)
        merge!(extra_metadata, eventlist.meta.extra)
    end

    return LightCurveMetadata(
        telescope,
        instrument,
        object,
        Float64(mjdref),
        (Float64(start_time), Float64(stop_time)),
        Float64(binsize),
        headers,  # Preserve ALL original headers
        extra_metadata,
    )
end
"""
    create_lightcurve(
        eventlist::EventList{TimeType, MetaType}, 
        binsize::Real;
        err_method::Symbol=:poisson,
        gaussian_errors::Union{Nothing,Vector{<:Real}}=nothing,
        tstart::Union{Nothing,Real}=nothing,
        tstop::Union{Nothing,Real}=nothing,
        energy_filter::Union{Nothing,Tuple{Real,Real}}=nothing,
        event_filter::Union{Nothing,Function}=nothing
    ) where {TimeType<:AbstractVector, MetaType<:FITSMetadata}

Create a light curve from an event list with comprehensive filtering and error handling.

This is the main function for creating light curves from X-ray event data. It supports 
comprehensive filtering options, multiple error calculation methods, and produces 
fully-documented light curve structures with complete metadata preservation.

# Arguments
- `eventlist::EventList{TimeType, MetaType}`: The input event list from `readevents`
- `binsize::Real`: Time bin size in seconds (must be positive)

# Keyword Arguments
- `err_method::Symbol=:poisson`: Error calculation method (`:poisson` or `:gaussian`)
- `gaussian_errors::Union{Nothing,Vector{<:Real}}=nothing`: User-provided errors (required for `:gaussian`)
- `tstart::Union{Nothing,Real}=nothing`: Start time for filtering (or `nothing` for data minimum)
- `tstop::Union{Nothing,Real}=nothing`: Stop time for filtering (or `nothing` for data maximum)
- `energy_filter::Union{Nothing,Tuple{Real,Real}}=nothing`: Energy range as `(emin, emax)` tuple in keV
- `event_filter::Union{Nothing,Function}=nothing`: Custom filter function taking EventList, returning boolean mask

# Returns
`LightCurve{T}`: Complete light curve structure with:
- Time-binned photon counts and statistical uncertainties
- Comprehensive metadata including processing history
- Additional properties (e.g., mean energy per bin)
- Full preservation of original FITS headers

# Error Methods
- `:poisson`: Uses Poisson statistics (σ = √N, with σ = 1 for N = 0)
- `:gaussian`: Uses user-provided Gaussian errors (must provide `gaussian_errors`)

# Filtering Options
1. **Time filtering**: Applied via `tstart` and `tstop` parameters
2. **Energy filtering**: Applied via `energy_filter` tuple (inclusive lower, exclusive upper)
3. **Custom filtering**: Applied via `event_filter` function for complex selection criteria

# Examples
```julia
# Basic usage with 1-second bins
ev = readevents("events.fits")
lc = create_lightcurve(ev, 1.0)
println("Created light curve with \$(length(lc)) bins")

# Energy-filtered light curve (0.5-10 keV)
lc_filtered = create_lightcurve(ev, 1.0, energy_filter=(0.5, 10.0))

# Time and energy filtering combined
lc_subset = create_lightcurve(ev, 1.0, 
                             tstart=1000.0, tstop=2000.0, 
                             energy_filter=(2.0, 8.0))

# Custom error calculation
expected_errs = sqrt.(expected_counts)  # Your theoretical errors
lc_custom = create_lightcurve(ev, 1.0, 
                             err_method=:gaussian, 
                             gaussian_errors=expected_errs)

# Complex custom filtering
function quality_filter(eventlist)
    # Example: filter based on multiple criteria
    return (eventlist.energies .> 0.3) .& 
           (eventlist.energies .< 12.0) .&
           (eventlist.pi .> 30)  # Assuming PI column exists
end

lc_quality = create_lightcurve(ev, 1.0, event_filter=quality_filter)

# High-resolution sub-second binning
lc_fast = create_lightcurve(ev, 0.1)  # 100ms bins
```

# Output Structure
The returned `LightCurve` provides:
- `lc.timebins`: Time bin centers
- `lc.counts`: Photon counts per bin
- `lc.count_error`: Statistical uncertainties
- `lc.exposure`: Exposure time per bin
- `lc.properties`: Additional derived properties (e.g., mean energy)
- `lc.metadata`: Complete observational and processing metadata

# Throws
- `ArgumentError`: If event list is empty
- `ArgumentError`: If bin size is not positive
- `ArgumentError`: If unsupported error method specified
- `ArgumentError`: If `:gaussian` method used without providing `gaussian_errors`
- `ArgumentError`: If `gaussian_errors` length doesn't match number of bins after filtering
- `ArgumentError`: If custom `event_filter` doesn't return boolean vector of correct length
- `ArgumentError`: If no events remain after any filtering step

# Performance Notes
- Uses vectorized operations for optimal performance with large event lists
- Memory-efficient binning algorithms from StatsBase.jl
- Filters applied in optimal order (energy first, then time) to minimize processing
- Type-stable implementation preserving input precision

# Implementation Details
The function performs these steps in order:
1. Input validation and type conversion
2. Custom event filtering (if specified)
3. Energy filtering (if specified)  
4. Time filtering (if specified)
5. Time bin creation with proper boundary handling
6. Event binning using optimized histogram algorithms
7. Error calculation based on specified method
8. Additional property calculation (mean energy, etc.)
9. Metadata extraction and preservation
10. Light curve structure creation

See also [`rebin`](@ref), [`LightCurve`](@ref), [`EventList`](@ref).
"""
function create_lightcurve(
    eventlist::EventList{TimeType,MetaType},
    binsize::Real;
    err_method::Symbol = :poisson,
    gaussian_errors::Union{Nothing,Vector{<:Real}} = nothing,
    tstart::Union{Nothing,Real} = nothing,
    tstop::Union{Nothing,Real} = nothing,
    energy_filter::Union{Nothing,Tuple{Real,Real}} = nothing,
    event_filter::Union{Nothing,Function} = nothing,
) where {TimeType<:AbstractVector,MetaType<:FITSMetadata}

    # Extract the element type from the vector type
    T = eltype(TimeType)

    # Validate all inputs first (but not gaussian_errors length yet)
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
    filtered_times, filtered_energies, start_time, stop_time =
        apply_event_filters(times, energies, tstart, tstop, energy_filter)

    # Create time bins
    bin_edges, bin_centers = create_time_bins(start_time, stop_time, binsize_t)

    # Bin the events
    counts = bin_events(filtered_times, bin_edges)

    # CRITICAL: Validate gaussian_errors length IMMEDIATELY after binning
    # This must happen BEFORE any success messages or further processing
    if err_method === :gaussian && !isnothing(gaussian_errors)
        if length(gaussian_errors) != length(counts)
            throw(
                ArgumentError(
                    "Length of gaussian_errors ($(length(gaussian_errors))) must match number of bins ($(length(counts)))",
                ),
            )
        end
    end

    @info "Created light curve: $(length(bin_centers)) bins, bin size = $(binsize_t) s"

    # Calculate exposures and errors
    exposure = fill(binsize_t, length(bin_centers))
    errors =
        calculate_errors(counts, err_method, exposure; gaussian_errors = gaussian_errors)

    # Calculate additional properties
    properties = calculate_additional_properties(
        filtered_times,
        filtered_energies,
        bin_edges,
        bin_centers,
    )

    # Extract metadata - Fixed to work with EventList structure
    metadata = extract_metadata(
        eventlist,
        start_time,
        stop_time,
        binsize_t,
        filtered_times,
        energy_filter,
    )

    return LightCurve{T}(
        bin_centers,
        bin_edges,
        counts,
        errors,
        exposure,
        properties,
        metadata,
        err_method,
    )
end

"""
    rebin(lc::LightCurve{T}, new_binsize::Real; 
          gaussian_errors::Union{Nothing,Vector{T}}=nothing) where T

Rebin a light curve to a new (larger) time resolution with proper error propagation.

This function combines adjacent time bins to create a light curve with lower time 
resolution. It properly handles count accumulation, error propagation, and 
property averaging while preserving all metadata and processing history.

# Arguments
- `lc::LightCurve{T}`: Input light curve to rebin
- `new_binsize::Real`: New (larger) bin size in seconds

# Keyword Arguments
- `gaussian_errors::Union{Nothing,Vector{T}}=nothing`: New error values for rebinned curve
  (required if original light curve used `:gaussian` error method)

# Returns
`LightCurve{T}`: Rebinned light curve with updated metadata

# Rebinning Process
1. **Count accumulation**: Counts from multiple old bins are summed into new bins
2. **Error propagation**: 
   - Poisson: σ² = Σ(σᵢ²) → σ = √(Σ counts)  
   - Gaussian: Must provide new errors via `gaussian_errors` parameter
3. **Property averaging**: Properties weighted by counts in original bins
4. **Metadata update**: Preserves original information, adds rebinning history

# Examples
```julia
# Create original 1-second light curve
ev = readevents("events.fits")
lc1 = create_lightcurve(ev, 1.0)
println("Original: \$(length(lc1)) bins of 1.0 s")

# Rebin to 10-second resolution
lc10 = rebin(lc1, 10.0)
println("Rebinned: \$(length(lc10)) bins of 10.0 s")

# Rebin with custom Gaussian errors
new_errors = sqrt.(expected_counts_10s)  # Your new error estimates
lc10_custom = rebin(lc1, 10.0, gaussian_errors=new_errors)

# Multiple rebinning steps
lc100 = rebin(lc10, 100.0)  # 1s → 10s → 100s
println("Final: \$(length(lc100)) bins of 100.0 s")

# Access rebinning history
println("Original bin size: ", lc100.metadata.extra["original_binsize"])
```

# Constraints
- `new_binsize` must be larger than current bin size
- For light curves with `:gaussian` errors, must provide `gaussian_errors`
- Properties are weighted-averaged (empty bins get zero values)
- Time alignment preserved from original binning

# Throws
- `ArgumentError`: If `new_binsize ≤ current_bin_size`
- `ArgumentError`: If `:gaussian` error method without providing `gaussian_errors`
- `ArgumentError`: If `gaussian_errors` length doesn't match number of new bins

# Performance Notes
- Efficient vectorized operations for large light curves
- Memory-efficient bin assignment using integer arithmetic
- Minimal memory allocation through pre-allocated arrays

# Statistical Considerations
- Rebinning reduces time resolution but improves signal-to-noise
- Count statistics remain valid (Poisson → Poisson)
- Properties may lose fine-scale variability information
- Metadata preserves full processing chain for reproducibility

See also [`create_lightcurve`](@ref), [`LightCurve`](@ref).
"""
function rebin(
    lc::LightCurve{T},
    new_binsize::Real;
    gaussian_errors::Union{Nothing,Vector{T}} = nothing,
) where {T}
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

    new_edges = [start_bin + i * new_binsize_t for i = 0:num_bins]
    new_centers = [start_bin + (i + 0.5) * new_binsize_t for i = 0:(num_bins-1)]

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
        throw(
            ArgumentError(
                "Gaussian errors must be provided when rebinning a light curve with Gaussian errors",
            ),
        )
    end

    new_errors = calculate_errors(
        new_counts,
        lc.err_method,
        new_exposure;
        gaussian_errors = gaussian_errors,
    )

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
            Dict{String,Any}("original_binsize" => Float64(old_binsize)),
        ),
    )

    return LightCurve{T}(
        new_centers,
        new_edges,
        new_counts,
        new_errors,
        new_exposure,
        new_properties,
        new_metadata,
        lc.err_method,
    )
end

# Array interface implementations with documentation
"""
    length(lc::LightCurve)

Return the number of time bins in the light curve.

# Examples
```julia
lc = create_lightcurve(ev, 1.0)
println("Light curve has \$(length(lc)) time bins")
```
"""
Base.length(lc::LightCurve) = length(lc.timebins)

"""
    size(lc::LightCurve)

Return the dimensions of the light curve as a tuple (for array interface compatibility).

# Examples
```julia
lc = create_lightcurve(ev, 1.0)
println("Light curve size: \$(size(lc))")  # (n_bins,)
```
"""
Base.size(lc::LightCurve) = (length(lc.timebins),)

"""
    getindex(lc::LightCurve, i::Int)

Get a (time, counts) tuple for the i-th time bin.

# Examples
```julia
lc = create_lightcurve(ev, 1.0)
time, counts = lc[1]  # First bin
println("Bin 1: time=\$time, counts=\$counts")
```
"""
Base.getindex(lc::LightCurve, i::Int) = (lc.timebins[i], lc.counts[i])

"""
    getindex(lc::LightCurve, r::UnitRange{Int})

Get (time, counts) tuples for a range of time bins.

# Examples
```julia
lc = create_lightcurve(ev, 1.0)
first_five = lc[1:5]  # First 5 bins as vector of tuples
```
"""
Base.getindex(lc::LightCurve, r::UnitRange{Int}) =
    [(lc.timebins[i], lc.counts[i]) for i in r]

"""
    iterate(lc::LightCurve)

Enable iteration over light curve bins, yielding (time, counts) tuples.

# Examples
```julia
lc = create_lightcurve(ev, 1.0)
for (time, counts) in lc
    println("Time: \$time, Counts: \$counts")
end
```
"""
Base.iterate(lc::LightCurve) =
    isempty(lc.timebins) ? nothing : ((lc.timebins[1], lc.counts[1]), 2)

"""
    iterate(lc::LightCurve, state)

Continue iteration over light curve bins.
"""
Base.iterate(lc::LightCurve, state) =
    state > length(lc.timebins) ? nothing :
    ((lc.timebins[state], lc.counts[state]), state + 1)
