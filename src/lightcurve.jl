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

A parametric structure for storing computed event properties per time bin.

$TYPEDEF

# Description
`EventProperty{T}` holds additional event properties beyond basic time and energy measurements.
These properties are typically calculated quantities such as mean energy, hardness ratios, 
or other derived statistics computed for each time bin of an observation.

# Type Parameters
- `T`: The numeric type of the property values (e.g., `Float64`, `Int64`)

# Fields
$FIELDS

# Constructors
```julia
EventProperty{T}(name::Symbol, values::Vector{T}, unit::String) where T
```

# Examples
```julia
# Mean energy per time bin
mean_energy = EventProperty{Float64}(:mean_energy, [1.2, 1.5, 1.8], "keV")

# Photon count rates
count_rate = EventProperty{Float64}(:rate, [10.5, 12.1, 9.8], "counts/s")

# Hardness ratio (dimensionless)
hardness = EventProperty{Float64}(:hardness, [0.3, 0.45, 0.2], "")

# Access properties
println(mean_energy.name)     # :mean_energy
println(mean_energy.values)   # [1.2, 1.5, 1.8]
println(mean_energy.unit)     # "keV"
```

# See Also
- [`LightCurveMetadata`](@ref): Metadata structure for light curves
- [`create_lightcurve`](@ref): Function to create light curves with properties
"""
struct EventProperty{T}
    "Name identifier for the property"
    name::Symbol
    "Property values for each time bin"
    values::Vector{T}
    "Physical units of the property values (empty string for dimensionless)"
    unit::String
end

"""
    LightCurveMetadata

Comprehensive metadata container for astronomical light curves.

$TYPEDEF

# Description
`LightCurveMetadata` stores all essential information about light curve creation and processing,
including observational parameters, timing information, and data provenance. This structure
ensures complete traceability and reproducibility of light curve analysis.

# Fields
$FIELDS

# Constructors
```julia
LightCurveMetadata(telescope, instrument, object, mjdref, time_range, 
                   bin_size, headers, extra)
```

# Usage Notes
- `mjdref` serves as the time reference point for all relative time measurements
- `time_range` uses the same time system as the light curve data
- `headers` preserves original FITS metadata for data provenance
- `extra` allows storage of custom processing parameters and derived information

# Examples
```julia
# Typical metadata creation (usually done automatically)
metadata = LightCurveMetadata(
    "NICER",                           # telescope
    "XTI",                             # instrument  
    "Crab Pulsar",                     # object
    58000.0,                           # mjdref
    (1000.0, 5000.0),                 # time_range
    1.0,                              # bin_size (seconds)
    [Dict("TELESCOP" => "NICER")],    # headers
    Dict("filter" => "0.3-10 keV")   # extra
)

# Access metadata fields
println(metadata.telescope)     # "NICER"
println(metadata.bin_size)      # 1.0
println(metadata.time_range)    # (1000.0, 5000.0)
println(metadata.extra["filter"]) # "0.3-10 keV"

# Check observation duration
duration = metadata.time_range[2] - metadata.time_range[1]  # 4000.0 seconds
```

# Implementation Details
- All time values should be in the same units (typically seconds)
- `headers` vector allows multiple FITS extensions to be preserved
- `extra` dictionary is extensible for mission-specific or analysis-specific metadata

# See Also
- [`EventProperty`](@ref): Structure for storing computed event properties
- [`create_lightcurve`](@ref): Primary function for light curve creation
"""
struct LightCurveMetadata
    "Name of the telescope/mission (e.g., 'NICER', 'XMM-Newton')"
    telescope::String
    "Name of the instrument (e.g., 'XTI', 'EPIC-pn')"
    instrument::String
    "Name or identifier of the observed astronomical object"
    object::String
    "Modified Julian Date reference time (MJD)"
    mjdref::Float64
    "Start and stop times of the light curve (relative to mjdref)"
    time_range::Tuple{Float64,Float64}
    "Time bin size in seconds"
    bin_size::Float64
    "Original FITS headers preserving data provenance"
    headers::Vector{Dict{String,Any}}
    "Additional metadata and processing information"
    extra::Dict{String,Any}
end

"""
    LightCurve{T} <: AbstractLightCurve{T}

Main structure for binned time series with statistical uncertainties and derived properties.

$TYPEDEF

# Description
`LightCurve{T}` represents binned photon count data with associated uncertainties, exposure times,
and computed properties. This is the primary data structure for time series analysis in
X-ray astronomy.

# Type Parameters
- `T`: Numeric type for time and derived quantities (typically `Float64`)

# Fields
$FIELDS

# Interface
- `length(lc)`: Number of time bins
- `lc[i]`: Get (time, counts) tuple for bin i
- Array-like indexing and iteration support

# Examples
```julia
# Create from event list
ev = readevents("events.fits")
lc = create_lightcurve(ev, 1.0)  # 1-second bins

# Access data
println("Time bins: ", lc.time[1:5])
println("Counts: ", lc.counts[1:5])
println("Errors: ", lc.count_error[1:5])

# Basic operations
println("Total counts: ", sum(lc.counts))
println("Mean count rate: ", mean(lc.counts ./ lc.exposure))
```

# See Also
- [`create_lightcurve`](@ref): Primary constructor function
- [`rebin`](@ref): Time binning operations
- [`EventProperty`](@ref): Additional computed properties
"""
mutable struct LightCurve{T} <: AbstractLightCurve{T}
    "Time bin centers"
    time::Vector{T}
    "Time bin width (scalar) or individual bin widths (vector)"
    dt::Union{T,Vector{T}}
    "Photon counts in each bin"
    counts::Vector{Int}
    "Statistical uncertainties on counts"
    count_error::Union{Nothing,Vector{T}}
    "Exposure time for each bin"
    exposure::Union{Nothing,Vector{T}}
    "Additional computed properties"
    properties::Vector{EventProperty{T}}
    "Comprehensive metadata"
    metadata::LightCurveMetadata
    "Error calculation method (:poisson or :gaussian)"
    err_method::Symbol
end

"""
    calculate_errors(counts, method, exposure; gaussian_errors=nothing)

Calculate statistical uncertainties for count data using vectorized operations.

$SIGNATURES

# Arguments
- `counts::Vector{Int}`: Photon counts in each bin
- `method::Symbol`: Error calculation method (`:poisson` or `:gaussian`)
- `exposure::Vector{T}`: Exposure times (interface compatibility)

# Keyword Arguments
- `gaussian_errors::Union{Nothing,Vector{T}}`: User-provided errors for `:gaussian` method

# Returns
`Vector{T}`: Statistical uncertainties for each bin

# Methods
- `:poisson`: σ = √N (σ = 1 for N = 0)
- `:gaussian`: Uses provided gaussian_errors

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
```

# Throws
- `ArgumentError`: Invalid method or missing gaussian_errors
"""

function calculate_errors(
    counts::Vector{Int},
    method::Symbol;
    gaussian_errors::Union{Nothing,Vector{<:Real}} = nothing
)
    if method === :poisson
        return @. sqrt(max(counts, 1))
    elseif method === :gaussian
        isnothing(gaussian_errors) && throw(ArgumentError(
            "Gaussian errors must be provided when using :gaussian method"
        ))
        length(gaussian_errors) != length(counts) && throw(ArgumentError(
            "Length of gaussian_errors must match length of counts"
        ))
        return convert(Vector{Float64}, gaussian_errors)
    else
        throw(ArgumentError("Unsupported error method: $method. Use :poisson or :gaussian"))
    end
end
"""
    set_errors!(lc::LightCurve{T}) -> LightCurve{T}
    set_errors!(lc::LightCurve{T}, errors::Vector{<:Real}) -> LightCurve{T}
Set error calculation method and update uncertainties in-place.
$SIGNATURES
# Methods
- `set_errors!(lc)`: Sets Poisson errors (σ = √N)
- `set_errors!(lc, errors)`: Sets custom Gaussian errors
# Examples
```julia
# Set Poisson errors
set_errors!(lc)
# Set custom errors
custom_errors = [0.5, 1.2, 0.8, 1.1]
set_errors!(lc, custom_errors)
```
# Throws
- `ArgumentError`: If custom errors length doesn't match bin count
"""
function set_errors!(lc::LightCurve{T}) where T
    lc.err_method = :poisson
    lc.count_error = convert(Vector{T}, calculate_errors(lc.counts, :poisson))
    return lc
end

function set_errors!(lc::LightCurve{T}, errors::Vector{<:Real}) where T
    length(errors) != length(lc.counts) && throw(ArgumentError(
        "Length of errors must match number of bins"
    ))
    lc.err_method = :gaussian
    lc.count_error = convert(Vector{T}, errors)
    return lc
end

"""
    calculate_errors!(lc::LightCurve{T}) -> Vector{T}
Recalculate and update uncertainties based on current error method.
$SIGNATURES
Updates the `count_error` field of the light curve using the current `err_method` 
and `counts` data. This is useful when the counts have been modified but you want 
to maintain the same error calculation approach.
# Examples
```julia
# After modifying counts, recalculate errors using existing method
lc.counts .+= background_correction
calculate_errors!(lc)  # Recalculates using lc.err_method
```
# Returns
- `Vector{T}`: Updated error vector (also stored in `lc.count_error`)
# See Also
- [`set_errors!`](@ref): Set error method and calculate errors
"""
function calculate_errors!(lc::LightCurve{T}) where T
    errors = calculate_errors(lc.counts, lc.err_method)
    lc.count_error = convert(Vector{T}, errors)
    return lc.count_error
end
"""
    create_time_bins(start_time, stop_time, binsize) -> (bin_edges, bin_centers)

Create uniform time bin edges and centers for light curve binning.

$SIGNATURES

Creates a uniform time grid covering the specified range with proper alignment
to bin boundaries and complete data coverage.

# Arguments
- `start_time::T`: Start time of the data range
- `stop_time::T`: Stop time of the data range  
- `binsize::T`: Size of each time bin

# Returns
- `bin_edges`: Array of bin boundaries (length = n_bins + 1)
- `bin_centers`: Array of bin center times (length = n_bins)

# Examples
```julia
# 1-second bins from 1000 to 1100 seconds
edges, centers = create_time_bins(1000.0, 1100.0, 1.0)
println(length(edges))    # 101 (100 bins + 1)
println(centers[1])       # 1000.5

# Sub-second binning
edges, centers = create_time_bins(0.0, 10.0, 0.1)
println(length(centers))  # 100
```

# Implementation
- Aligns start_bin to bin_size boundaries for consistency
- Ensures complete coverage including stop_time
- Centers calculated as bin_start + 0.5 * bin_size
"""
function create_time_bins(start_time::T, stop_time::T, binsize::T) where T
    start_bin = floor(start_time / binsize) * binsize
    time_span = stop_time - start_bin
    num_bins = max(1, ceil(Int, time_span / binsize))
    
    if start_bin + num_bins * binsize <= stop_time
        num_bins += 1
    end
    
    edges = [start_bin + i * binsize for i = 0:num_bins]
    centers = [start_bin + (i + 0.5) * binsize for i = 0:(num_bins-1)]
    
    return edges, centers
end

"""
    bin_events(times, bin_edges) -> Vector{Int}

Bin event arrival times into histogram counts using optimized algorithms.

$SIGNATURES

Efficiently bins photon arrival times using StatsBase.jl with the rightmost
bin edge made inclusive to capture all events.

# Arguments
- `times::Vector{T}`: Event arrival times (need not be sorted)
- `bin_edges::Vector{T}`: Time bin edges (length = n_bins + 1)

# Returns
`Vector{Int}`: Photon counts in each bin

# Examples
```julia
times = [1.1, 1.3, 1.7, 2.2, 2.8, 3.1]
edges = [1.0, 2.0, 3.0, 4.0]  # 3 bins

counts = bin_events(times, edges)
# Result: [3, 2, 1]
```

# Binning Rules
- All bins except last are half-open [a, b)
- Last bin is closed [a, b] to capture boundary events
- Events on interior boundaries belong to the right bin

# Throws
- `ArgumentError`: If fewer than 2 bin edges provided
"""
function bin_events(times::AbstractVector, dt::Vector{T}) where T
    length(dt) < 2 && throw(ArgumentError("Need at least 2 bin edges"))
    
    adjusted_edges = copy(dt)
    if length(adjusted_edges) > 1
        adjusted_edges[end] = nextfloat(adjusted_edges[end])
    end
    
    hist = fit(Histogram, times, adjusted_edges)
    return Vector{Int}(hist.weights)
end

"""
    apply_filters(times, energies, tstart, tstop, energy_filter) -> (filtered_times, filtered_energies, start_t, stop_t)

Apply time and energy filters to event data with validation.

$SIGNATURES

Applies filtering operations in optimal order with comprehensive validation
and automatic time range determination.

# Arguments
- `times::Vector{T}`: Event arrival times
- `energies::Union{Nothing,Vector{T}}`: Event energies (optional)
- `tstart::Union{Nothing,Real}`: Start time filter (optional)
- `tstop::Union{Nothing,Real}`: Stop time filter (optional)
- `energy_filter::Union{Nothing,Tuple{Real,Real}}`: Energy range (emin, emax)

# Returns
- Filtered times and energies
- Final start and stop times used

# Examples
```julia
# Energy filter only
filtered_times, filtered_energies, start_t, stop_t = apply_filters(
    times, energies, nothing, nothing, (0.5, 10.0)
)

# Both time and energy filters
filtered_times, filtered_energies, start_t, stop_t = apply_filters(
    times, energies, 1000.0, 2000.0, (0.5, 10.0)
)
```

# Throws
- `ArgumentError`: If no events remain after filtering
"""
function apply_filters(
    times::AbstractVector{T},
    energies::Union{Nothing,AbstractVector{T}},
    tstart::Union{Nothing,Real},
    tstop::Union{Nothing,Real},
    energy_filter::Union{Nothing,Tuple{Real,Real}}
) where T
    # Start with all indices
    mask = trues(length(times))
    
    # Apply energy filter first if provided and energies exist
    if !isnothing(energy_filter) && !isnothing(energies)
        emin, emax = energy_filter
        energy_mask = (energies .>= emin) .& (energies .< emax)
        mask = mask .& energy_mask
    end
    
    # Apply time filters
    if !isnothing(tstart)
        mask = mask .& (times .>= tstart)
    end
    if !isnothing(tstop)
        mask = mask .& (times .<= tstop)
    end
    
    # Check if any events remain
    if !any(mask)
        throw(ArgumentError("No events remain after applying filters"))
    end
    
    # Apply the mask
    filtered_times = times[mask]
    filtered_energies = isnothing(energies) ? nothing : energies[mask]
    
    # Calculate time range
    start_t = isnothing(tstart) ? minimum(times) : tstart
    stop_t = isnothing(tstop) ? maximum(times) : tstop
    
    return filtered_times, filtered_energies, start_t, stop_t
end
"""
    apply_filters(times, energies, tstart, tstop, energy_filter)

Basic event filtering without GTI consideration.

# Arguments
- `times::AbstractVector{T}`: Event arrival times
- `energies::Union{Nothing,AbstractVector{T}}`: Event energies (optional)
- `tstart::Union{Nothing,Real}`: Start time filter (inclusive)
- `tstop::Union{Nothing,Real}`: Stop time filter (inclusive)
- `energy_filter::Union{Nothing,Tuple{Real,Real}}`: Energy range (emin, emax)

# Returns
`Tuple{Vector, Union{Nothing,Vector}, Real, Real}`: 
- Filtered times
- Filtered energies (if provided)
- Actual start time
- Actual stop time

# Examples
```julia
# Time filtering only
filtered_times, _, start_t, stop_t = apply_filters(times, nothing, 1000.0, 2000.0, nothing)

# Energy and time filtering
filtered_times, filtered_energies, start_t, stop_t = apply_filters(
    times, energies, 1000.0, 2000.0, (0.5, 10.0)
)
Notes
- Energy filter is applied as: emin ≤ energy < emax
- Time filter is applied as: tstart ≤ time ≤ tstop
"""
function apply_filters(
    times::AbstractVector{T},
    energies::Union{Nothing,AbstractVector{T}},
    eventlist::EventList,
    tstart::Union{Nothing,Real},
    tstop::Union{Nothing,Real},
    energy_filter::Union{Nothing,Tuple{Real,Real}},
    binsize::Real
) where T
    mask = trues(length(times))
    
    # Apply energy filter
    if !isnothing(energy_filter) && !isnothing(energies)
        emin, emax = energy_filter
        mask = mask .& (energies .>= emin) .& (energies .< emax)
    end
    
    # Apply time filters
    if !isnothing(tstart)
        mask = mask .& (times .>= tstart)
    end
    if !isnothing(tstop)
        mask = mask .& (times .<= tstop)
    end
    # If GTI is present, apply GTI mask
    if has_gti(eventlist)
        gti_mask, _ = create_gti_mask(times, eventlist.meta.gti, dt=binsize)
        mask = mask .& gti_mask
    end
    
    !any(mask) && throw(ArgumentError("No events remain after applying filters"))
    
    filtered_times = times[mask]
    filtered_energies = isnothing(energies) ? nothing : energies[mask]
    start_t = isnothing(tstart) ? minimum(filtered_times) : tstart
    stop_t = isnothing(tstop) ? maximum(filtered_times) : tstop
    
    return filtered_times, filtered_energies, start_t, stop_t
end

"""
    calculate_event_properties(times, energies, dt, bin_centers) -> Vector{EventProperty}

Calculate derived properties for each time bin from event data.

$SIGNATURES

Computes additional properties like mean energy per bin from the input events.
Currently calculates mean energy when energy data is available.

# Arguments
- `times::AbstractVector`: Event arrival times
- `energies::Union{Nothing,AbstractVector}`: Event energies (optional)
- `dt::Vector{T}`: Time bin edges
- `bin_centers::Vector{T}`: Time bin centers

# Returns
`Vector{EventProperty}`: Array of computed properties

# Examples
```julia
properties = calculate_event_properties(times, energies, edges, centers)
# Returns EventProperty with mean_energy if energies provided
```

# Properties Calculated
- `:mean_energy`: Average energy per bin (when energies available)
"""
function calculate_event_properties(
    times::AbstractVector,
    energies::Union{Nothing,AbstractVector},
    dt::Vector{T},
    bin_centers::Vector{T}
) where T
    properties = Vector{EventProperty}()
    
    if !isnothing(energies) && !isempty(energies) && length(bin_centers) > 0
        start_bin = dt[1]
        binsize = length(bin_centers) == 1 ? 
            (length(dt) > 1 ? dt[2] - dt[1] : T(1)) :
            bin_centers[2] - bin_centers[1]
        
        energy_sums = zeros(T, length(bin_centers))
        energy_counts = zeros(Int, length(bin_centers))
        
        for (t, e) in zip(times, energies)
            bin_idx = floor(Int, (t - start_bin) / binsize) + 1
            if 1 <= bin_idx <= length(bin_centers)
                energy_sums[bin_idx] += T(e)
                energy_counts[bin_idx] += 1
            end
        end
        
        mean_energy = @. ifelse(energy_counts > 0, energy_sums / energy_counts, zero(T))
        push!(properties, EventProperty{T}(:mean_energy, mean_energy, "keV"))
    end
    
    return properties
end
"""
    extract_metadata(eventlist, start_time, stop_time, binsize, n_filtered_events, energy_filter) -> LightCurveMetadata

Extract and construct metadata from an EventList for light curve creation.

$SIGNATURES

Extracts observational metadata from FITS headers and constructs a comprehensive
metadata structure with filtering and processing information.

# Arguments
- `eventlist::EventList`: Source event list
- `start_time, stop_time`: Time range of the light curve
- `binsize`: Time bin size used
- `n_filtered_events`: Number of events after filtering
- `energy_filter`: Energy filter applied (if any)

# Returns
`LightCurveMetadata`: Complete metadata structure

# Implementation
- Handles various header formats (FITSHeader, Vector, Dict)
- Extracts standard FITS keywords (TELESCOP, INSTRUME, OBJECT, MJDREF)
- Preserves filtering and processing history
"""
function extract_metadata(eventlist::EventList, start_time, stop_time, binsize, n_filtered_events, energy_filter)
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
    
    telescope = get(first_header, "TELESCOP", get(first_header, "TELESCOPE", ""))
    instrument = get(first_header, "INSTRUME", get(first_header, "INSTRUMENT", ""))
    object = get(first_header, "OBJECT", get(first_header, "TARGET", ""))
    mjdref = get(first_header, "MJDREF", get(first_header, "MJDREFI", 0.0))
    
    # Handle both array of events and count of events
    n_events = if n_filtered_events isa AbstractVector
        length(n_filtered_events)
    else
        n_filtered_events  # Assume it's already a count
    end
    
    # Create extra metadata with GTI information[storing purpose]
    extra_metadata = Dict{String,Any}(
        "filtered_nevents" => n_events,
        "total_nevents" => length(eventlist.times),
        "energy_filter" => energy_filter,
        "binning_method" => "histogram",
        "gti" => eventlist.meta.gti  #GTI information[this 'eventlist.meta.gti' need to bw rembered since it is the one where u will all gti information:)]
    )
    
    if hasfield(typeof(eventlist.meta), :extra)
        merge!(extra_metadata, eventlist.meta.extra)
    end
    
    # Add GTI source information if available
    if !isnothing(eventlist.meta.gti_source)
        extra_metadata["gti_source"] = eventlist.meta.gti_source
    end
    
    return LightCurveMetadata(
        telescope, instrument, object, Float64(mjdref),
        (Float64(start_time), Float64(stop_time)), Float64(binsize),
        headers, extra_metadata
    )
end
"""
    create_lightcurve(eventlist, binsize; kwargs...) -> LightCurve

Create a binned light curve from an event list with filtering and error calculation.

$SIGNATURES

Main function for converting event lists to light curves with comprehensive
filtering options and automatic error calculation.

# Arguments
- `eventlist::EventList`: Input event data (should be pre-filtered for custom filters)
- `binsize::Real`: Time bin size in seconds

# Keyword Arguments
- `err_method::Symbol = :poisson`: Error calculation method (`:poisson` or `:gaussian`)
- `tstart::Union{Nothing,Real} = nothing`: Start time filter
- `tstop::Union{Nothing,Real} = nothing`: Stop time filter  
- `energy_filter::Union{Nothing,Tuple{Real,Real}} = nothing`: Energy range (emin, emax)

# Returns
`LightCurve{T}`: Binned light curve with metadata and computed properties

# Examples
```julia
# Basic light curve creation
lc = create_lightcurve(eventlist, 1.0)

# With time and energy filtering
lc = create_lightcurve(eventlist, 0.1; 
                      tstart=1000.0, tstop=2000.0,
                      energy_filter=(0.5, 10.0))

# With custom event filter (apply before lightcurve creation)
filtered_eventlist = filter_time(t -> t > 1000.0, eventlist)
lc = create_lightcurve(filtered_eventlist, 1.0)

# Multiple custom filters can be chained
filtered_eventlist = eventlist |>
    ev -> filter_time(t -> t > 1000.0, ev) |>
    ev -> filter_energy(e -> e < 10.0, ev)
lc = create_lightcurve(filtered_eventlist, 1.0)
```

# Notes
- Custom event filtering should be performed on the `EventList` before calling this function
- Use `filter_time!`, `filter_energy!`, or other filtering functions for custom filtering
- The function applies only basic time and energy filters specified in keyword arguments

# Throws
- `ArgumentError`: If event list is empty, binsize ≤ 0, invalid error method, or no events after filtering
"""
function create_lightcurve(
    eventlist::EventList{TimeType,MetaType},
    binsize::Real;
    err_method::Symbol = :poisson,
    tstart::Union{Nothing,Real} = nothing,
    tstop::Union{Nothing,Real} = nothing,
    energy_filter::Union{Nothing,Tuple{Real,Real}} = nothing
) where {TimeType<:AbstractVector, MetaType<:FITSMetadata}
    
    T = eltype(TimeType)
    
    # Validate inputs
    isempty(eventlist.times) && throw(ArgumentError("Event list is empty"))
    binsize <= 0 && throw(ArgumentError("Bin size must be positive"))
    !(err_method in [:poisson, :gaussian]) && throw(ArgumentError(
        "Unsupported error method: $err_method. Use :poisson or :gaussian"
    ))
    binsize_t = convert(T, binsize)
    

    filtered_times, filtered_energies, start_t, stop_t = apply_filters(
        eventlist.times,
        eventlist.energies,
        eventlist,
        tstart,
        tstop,
        energy_filter,
        binsize_t
    )
    
    isempty(filtered_times) && throw(ArgumentError("No events remain after filtering"))
    
    # Determine time range
    start_time = minimum(filtered_times)
    stop_time = maximum(filtered_times)
    
    # Create time bins and bin events
    dt, bin_centers = create_time_bins(start_time, stop_time, binsize_t)
    counts = bin_events(filtered_times, dt)
    
    @debug "Created light curve: $(length(bin_centers)) bins, bin size = $(binsize_t) s"
    
    # Calculate exposure and properties
    exposure = fill(binsize_t, length(bin_centers))
    properties = calculate_event_properties(filtered_times, filtered_energies, dt, bin_centers)
    
    # Extract metadata with GTI information
    actual_start = !isnothing(tstart) ? T(tstart) : start_time
    actual_stop = !isnothing(tstop) ? T(tstop) : stop_time
    metadata = extract_metadata(eventlist, actual_start, actual_stop, binsize_t, 
                              length(filtered_times), energy_filter)
    
    # Create light curve (errors will be calculated when needed)
    lc = LightCurve{T}(
        bin_centers, binsize_t, counts, nothing, exposure,
        properties, metadata, err_method
    )
    
    # Calculate initial errors
    calculate_errors!(lc)
    
    # Add debug info about GTI
    if has_gti(eventlist)
        @debug "GTI information preserved" n_intervals=size(eventlist.meta.gti, 1) time_range=extrema(eventlist.meta.gti)
    end
    
    return lc
end
"""
    rebin(lc::LightCurve, new_binsize::Real) -> LightCurve

Rebin a light curve to larger time bins while preserving total counts and properties.

$SIGNATURES

Creates a new light curve with larger time bins by combining adjacent bins.
Preserves total photon counts and recalculates weighted averages for properties.

# Arguments
- `lc::LightCurve`: Input light curve
- `new_binsize::Real`: New bin size (must be larger than current)

# Returns
`LightCurve{T}`: Rebinned light curve with updated metadata

# Examples
```julia
# Rebin from 0.1s to 1.0s bins
lc_1s = rebin(lc_0p1s, 1.0)

# Rebin to 10-second bins
lc_10s = rebin(lc, 10.0)
```

# Implementation
- Combines counts additively across bins
- Recalculates weighted averages for properties
- Updates metadata with original bin size information
- Preserves error calculation method

# Throws
- `ArgumentError`: If new bin size is not larger than current bin size
"""
function rebin(lc::LightCurve{T}, new_binsize::Real) where T
    new_binsize <= lc.metadata.bin_size && throw(ArgumentError(
        "New bin size must be larger than current bin size"
    ))
    
    old_binsize = T(lc.metadata.bin_size)
    new_binsize_t = convert(T, new_binsize)
    
    start_time = T(lc.metadata.time_range[1])
    stop_time = T(lc.metadata.time_range[2])
    
    # Create new bins
    new_edges, new_centers = create_time_bins(start_time, stop_time, new_binsize_t)
    
    # Rebin counts
    new_counts = zeros(Int, length(new_centers))
    start_bin = new_edges[1]
    
    for (i, time) in enumerate(lc.time)
        if lc.counts[i] > 0
            bin_idx = floor(Int, (time - start_bin) / new_binsize_t) + 1
            if 1 <= bin_idx <= length(new_counts)
                new_counts[bin_idx] += lc.counts[i]
            end
        end
    end
    
    # Handle properties
    new_properties = Vector{EventProperty}()
    for prop in lc.properties
        new_values = zeros(T, length(new_centers))
        counts = zeros(Int, length(new_centers))
        
        for (i, val) in enumerate(prop.values)
            if lc.counts[i] > 0
                bin_idx = floor(Int, (lc.time[i] - start_bin) / new_binsize_t) + 1
                if 1 <= bin_idx <= length(new_values)
                    new_values[bin_idx] += val * lc.counts[i]
                    counts[bin_idx] += lc.counts[i]
                end
            end
        end
        
        new_values = @. ifelse(counts > 0, new_values / counts, zero(T))
        push!(new_properties, EventProperty(prop.name, new_values, prop.unit))
    end
    
    # Update metadata
    new_metadata = LightCurveMetadata(
        lc.metadata.telescope, lc.metadata.instrument, lc.metadata.object,
        lc.metadata.mjdref, lc.metadata.time_range, Float64(new_binsize_t),
        lc.metadata.headers,
        merge(lc.metadata.extra, Dict{String,Any}("original_binsize" => Float64(old_binsize)))
    )
    
    # Create rebinned light curve
    rebinned_lc = LightCurve{T}(
        new_centers, new_binsize_t, new_counts, nothing,
        fill(new_binsize_t, length(new_centers)), new_properties,
        new_metadata, lc.err_method
    )
    
    # Calculate errors for rebinned curve
    calculate_errors!(rebinned_lc)
    
    return rebinned_lc
end
"""
Array interface implementation for LightCurve structures.

Provides standard Julia array operations for convenient access to time-count pairs.

# Interface Methods
- `length(lc)`: Number of time bins
- `size(lc)`: Tuple with light curve dimensions  
- `lc[i]`: Get (time, counts) tuple for bin i
- `lc[range]`: Get array of (time, counts) tuples for range
- `iterate(lc)`: Iterator support for loops

# Examples
```julia
# Basic access
println(length(lc))        # Number of bins
println(lc[1])            # (time, counts) for first bin
println(lc[1:5])          # First 5 bins

# Iteration
for (time, counts) in lc
    println("Time: $time, Counts: $counts")
end
```
"""
Base.length(lc::LightCurve) = length(lc.time)
Base.size(lc::LightCurve) = (length(lc.time),)
Base.lastindex(lc::LightCurve) = length(lc.time)
Base.getindex(lc::LightCurve, i::Int) = (lc.time[i], lc.counts[i])
Base.getindex(lc::LightCurve, r::UnitRange{Int}) = [(lc.time[i], lc.counts[i]) for i in r]

Base.iterate(lc::LightCurve) = 
    isempty(lc.time) ? nothing : ((lc.time[1], lc.counts[1]), 2)
Base.iterate(lc::LightCurve, state) = 
state > length(lc.time) ? nothing : ((lc.time[state], lc.counts[state]), state + 1)