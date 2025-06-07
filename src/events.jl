"""
    FITSMetadata{H}

Metadata associated with a FITS or events file.

# Fields
- `filepath::String`: Path to the FITS file
- `hdu::Int`: HDU index that the metadata was read from
- `energy_units::Union{Nothing,String}`: Units of energy (column name: ENERGY, PI, or PHA)
- `extra_columns::Dict{String,Vector}`: Extra columns that were requested during read
- `headers::H`: FITS headers from the selected HDU

# Examples
```julia
# Metadata is typically created automatically when reading events
ev = readevents("data.fits")
println(ev.meta.filepath)  # Shows the file path
println(ev.meta.energy_units)  # Shows "PI", "ENERGY", or "PHA"
```
"""
struct FITSMetadata{H}
    "Path to the FITS file"
    filepath::String
    "HDU index that the metadata was read from"
    hdu::Int
    "Units of energy (column name: ENERGY, PI, or PHA)"
    energy_units::Union{Nothing,String}
    "Extra columns that were requested during read"
    extra_columns::Dict{String,Vector}
    "FITS headers from the selected HDU"
    headers::H
end

function Base.show(io::IO, ::MIME"text/plain", m::FITSMetadata)
    println(
        io,
        "FITSMetadata for $(basename(m.filepath))[$(m.hdu)] with $(length(m.extra_columns)) extra column(s)",
    )
end

"""
    EventList{TimeType, MetaType <: FITSMetadata}

Container for an events list storing times, energies, and associated metadata.

# Fields
- `times::TimeType`: Vector with recorded times
- `energies::Union{Nothing,TimeType}`: Vector with recorded energies (or `nothing` if no energy data)
- `meta::MetaType`: Metadata from FITS file

# Constructors
```julia
# Read from FITS file (recommended)
ev = readevents("events.fits")

# Create directly for testing (simplified constructor)
ev = EventList([1.0, 2.0, 3.0], [0.5, 1.2, 2.1])  # times and energies
ev = EventList([1.0, 2.0, 3.0])  # times only
```

# Interface
- `length(ev)`: Number of events
- `times(ev)`: Access times vector
- `energies(ev)`: Access energies vector (may be `nothing`)
- `has_energies(ev)`: Check if energies are present

# Filtering
EventList supports composable filtering operations:
```julia
# Filter by time (in-place)
filter_time!(t -> t > 100.0, ev)

# Filter by energy (in-place)  
filter_energy!(energy_val -> energy_val < 10.0, ev)

# Non-mutating versions
ev_filtered = filter_time(t -> t > 100.0, ev)
ev_filtered = filter_energy(energy_val -> energy_val < 10.0, ev)

# Composable filtering
filter_energy!(x -> x < 10.0, filter_time!(t -> t > 100.0, ev))
```

Generally should not be directly constructed, but read from file using [`readevents`](@ref).
"""
struct EventList{TimeType<:AbstractVector,MetaType<:FITSMetadata}
    "Vector with recorded times"
    times::TimeType
    "Vector with recorded energies (else `nothing`)"
    energies::Union{Nothing,TimeType}
    "Metadata from FITS file"
    meta::MetaType
end

"""
    EventList(times::Vector{T}, energies::Union{Nothing,Vector{T}}=nothing) where T

Simple constructor for testing without FITS files. Creates an EventList with dummy metadata.

# Arguments
- `times::Vector{T}`: Vector of event times
- `energies::Union{Nothing,Vector{T}}`: Optional vector of event energies

# Examples
```julia
# Times only
ev = EventList([1.0, 2.0, 3.0])

# Times and energies
ev = EventList([1.0, 2.0, 3.0], [0.5, 1.2, 2.1])
```
"""
function EventList(times::Vector{T}, energies::Union{Nothing,Vector{T}} = nothing) where {T}
    dummy_meta = FITSMetadata(
        "",  # filepath
        1,   # hdu
        nothing,  # energy_units
        Dict{String,Vector}(),  # extra_columns
        Dict{String,Any}(),  # headers
    )
    EventList(times, energies, dummy_meta)
end

function Base.show(io::IO, ::MIME"text/plain", ev::EventList)
    print(io, "EventList with $(length(ev.times)) times")
    if !isnothing(ev.energies)
        print(io, " and energies")
    end
    println(io)
end

# ============================================================================
# Interface Methods
# ============================================================================

"""
    length(ev::EventList)

Return the number of events in the EventList.
"""
Base.length(ev::EventList) = length(ev.times)

"""
    size(ev::EventList)

Return the size of the EventList as a tuple.
"""
Base.size(ev::EventList) = (length(ev),)

"""
    times(ev::EventList)

Access the times vector of an EventList.

# Examples
```julia
ev = readevents("data.fits")
time_data = times(ev)  # Get times vector
```
"""
times(ev::EventList) = ev.times

"""
    energies(ev::EventList)

Access the energies vector of an EventList. Returns `nothing` if no energies are present.

# Examples
```julia
ev = readevents("data.fits")
energy_data = energies(ev)  # May be nothing
if !isnothing(energy_data)
    println("Energy range: \$(extrema(energy_data))")
end
```
"""
energies(ev::EventList) = ev.energies

"""
    has_energies(ev::EventList)

Check whether the EventList contains energy information.

# Examples
```julia
ev = readevents("data.fits")
if has_energies(ev)
    println("Energy data available")
end
```
"""
has_energies(ev::EventList) = !isnothing(ev.energies)

# ============================================================================
# Filtering Functions (Composable and In-Place)
# ============================================================================

"""
    filter_time!(f, ev::EventList)

Filter all columns of the EventList based on a predicate `f` applied to the times. 
Modifies the EventList in-place for efficiency.

# Arguments
- `f`: Predicate function that takes a time value and returns a Boolean
- `ev::EventList`: EventList to filter (modified in-place)

# Returns
The modified EventList (for chaining operations)

# Examples
```julia
# Filter only positive times
filter_time!(t -> t > 0, ev)

# Filter times greater than some minimum using function composition
min_time = 100.0
filter_time!(x -> x > min_time, ev)

# Chaining filters
filter_energy!(x -> x < 10.0, filter_time!(t -> t > 100.0, ev))
```

See also [`filter_energy!`](@ref), [`filter_time`](@ref).
"""
filter_time!(f, ev::EventList) = filter_on!(f, ev.times, ev)

"""
    filter_energy!(f, ev::EventList)

Filter all columns of the EventList based on a predicate `f` applied to the energies. 
Modifies the EventList in-place for efficiency.

# Arguments
- `f`: Predicate function that takes an energy value and returns a Boolean
- `ev::EventList`: EventList to filter (modified in-place)

# Returns
The modified EventList (for chaining operations)

# Examples
```julia
# Filter energies less than 10 keV
filter_energy!(energy_val -> energy_val < 10.0, ev)

# With function composition
max_energy = 10.0
filter_energy!(x -> x < max_energy, ev)

# Chaining with time filter
filter_energy!(x -> x < 10.0, filter_time!(t -> t > 100.0, ev))
```

# Throws
- `AssertionError`: If the EventList has no energy data

See also [`filter_time!`](@ref), [`filter_energy`](@ref).
"""
function filter_energy!(f, ev::EventList)
    @assert !isnothing(ev.energies) "No energies present in the EventList."
    filter_on!(f, ev.energies, ev)
end

"""
    filter_on!(f, src_col::AbstractVector, ev::EventList)

Internal function to filter EventList based on predicate applied to source column.
Uses efficient in-place filtering adapted from Base.filter! implementation.

This function maintains consistency across all columns (times, energies, extra_columns)
by applying the same filtering mask derived from the source column.

# Arguments
- `f`: Predicate function applied to elements of `src_col`
- `src_col::AbstractVector`: Source column to generate filtering mask from
- `ev::EventList`: EventList to filter

# Implementation Notes
- Uses `eachindex` for portable iteration over array indices
- Modifies arrays in-place using a two-pointer technique
- Resizes all arrays to final filtered length
- Maintains type stability with `::Bool` annotation on predicate result
"""
function filter_on!(f, src_col::AbstractVector, ev::EventList)
    @assert size(src_col) == size(ev.times) "Source column size must match times size"

    # Modified from Base.filter! implementation for multiple arrays
    # Use two pointers: i for reading, j for writing
    j = firstindex(ev.times)

    for i in eachindex(ev.times)
        predicate = f(src_col[i])::Bool

        if predicate
            # Copy elements to new position
            ev.times[j] = ev.times[i]

            if !isnothing(ev.energies)
                ev.energies[j] = ev.energies[i]
            end

            # Handle extra columns
            for (_, col) in ev.meta.extra_columns
                col[j] = col[i]
            end

            j = nextind(ev.times, j)
        end
    end

    # Resize all arrays to new length
    if j <= lastindex(ev.times)
        new_length = j - 1
        resize!(ev.times, new_length)

        if !isnothing(ev.energies)
            resize!(ev.energies, new_length)
        end

        for (_, col) in ev.meta.extra_columns
            resize!(col, new_length)
        end
    end

    ev
end

# ============================================================================
# Non-mutating Filter Functions
# ============================================================================

"""
    filter_time(f, ev::EventList)

Return a new EventList with events filtered by predicate `f` applied to times.
This is the non-mutating version of [`filter_time!`](@ref).

# Arguments
- `f`: Predicate function that takes a time value and returns a Boolean
- `ev::EventList`: EventList to filter (not modified)

# Returns
New EventList with filtered events

# Examples
```julia
# Create filtered copy
ev_filtered = filter_time(t -> t > 100.0, ev)

# Original EventList is unchanged
println(length(ev))  # Original length
println(length(ev_filtered))  # Filtered length
```

See also [`filter_time!`](@ref), [`filter_energy`](@ref).
"""
function filter_time(f, ev::EventList)
    new_ev = deepcopy(ev)
    filter_time!(f, new_ev)
end

"""
    filter_energy(f, ev::EventList)

Return a new EventList with events filtered by predicate `f` applied to energies.
This is the non-mutating version of [`filter_energy!`](@ref).

# Arguments
- `f`: Predicate function that takes an energy value and returns a Boolean
- `ev::EventList`: EventList to filter (not modified)

# Returns
New EventList with filtered events

# Examples
```julia
# Create filtered copy
ev_filtered = filter_energy(energy_val -> energy_val < 10.0, ev)

# Original EventList is unchanged
println(length(ev))  # Original length
println(length(ev_filtered))  # Filtered length
```

# Throws
- `AssertionError`: If the EventList has no energy data

See also [`filter_energy!`](@ref), [`filter_time`](@ref).
"""
function filter_energy(f, ev::EventList)
    new_ev = deepcopy(ev)
    filter_energy!(f, new_ev)
end

# ============================================================================
# File Reading Functions
# ============================================================================

"""
    colnames(file::AbstractString; hdu = 2)

Return a vector of all column names of a FITS file, reading from the specified HDU.

# Arguments
- `file::AbstractString`: Path to FITS file
- `hdu::Int`: HDU index to read from (default: 2, typical for event data)

# Returns
Vector of column names as strings

# Examples
```julia
cols = colnames("events.fits")
println(cols)  # ["TIME", "PI", "RAWX", "RAWY", ...]

# Check if energy column exists
if "ENERGY" in colnames("events.fits")
    println("Energy data available")
end
```
"""
function colnames(file::AbstractString; hdu = 2)
    FITS(file) do f
        selected_hdu = f[hdu]
        FITSIO.colnames(selected_hdu)
    end
end

"""
    read_energy_column(hdu; energy_alternatives = ["ENERGY", "PI", "PHA"], T = Float64)

Attempt to read the energy column of an HDU from a list of alternative names.

This function provides a robust way to read energy data from FITS files, as different
missions and instruments use different column names for energy information.

# Arguments
- `hdu`: FITS HDU object to read from
- `energy_alternatives::Vector{String}`: List of column names to try (default: ["ENERGY", "PI", "PHA"])
- `T::Type`: Type to convert energy data to (default: Float64)

# Returns
`(column_name, data)` tuple where:
- `column_name::Union{Nothing,String}`: Name of the column that was successfully read, or `nothing`
- `data::Union{Nothing,Vector{T}}`: Energy data as Vector{T}, or `nothing` if no column found

# Examples
```julia
FITS("events.fits") do f
    hdu = f[2]
    col_name, energy_data = read_energy_column(hdu)
    if !isnothing(energy_data)
        println("Found energy data in column: \$col_name")
        println("Energy range: \$(extrema(energy_data))")
    end
end
```

# Implementation Notes
- Tries columns in order until one is successfully read
- Uses case-insensitive matching for column names
- Handles read errors gracefully by trying the next column
- Separated from main reading function for testability and clarity
- Type-stable with explicit return type annotation
- Added case_sensitive=false parameter: This tells FITSIO.jl to use the old behavior for backward compatibility
"""
function read_energy_column(
    hdu;
    energy_alternatives::Vector{String} = ["ENERGY", "PI", "PHA"],
    T::Type = Float64,
)::Tuple{Union{Nothing,String},Union{Nothing,Vector{T}}}

    # Get actual column names from the file
    all_cols = FITSIO.colnames(hdu)
    
    for col_name in energy_alternatives
        # Find matching column name (case-insensitive)
        actual_col = findfirst(col -> uppercase(col) == uppercase(col_name), all_cols)
        
        if !isnothing(actual_col)
            actual_col_name = all_cols[actual_col]
            try
                # Use the actual column name from the file
                data = read(hdu, actual_col_name, case_sensitive=false)
                return actual_col_name, convert(Vector{T}, data)
            catch
                # If this column exists but can't be read, try the next one
                continue
            end
        end
    end

    return nothing, nothing
end


"""
    readevents(path; kwargs...)

Read an [`EventList`](@ref) from a FITS file. Will attempt to read an energy
column if one exists.

This is the primary function for loading X-ray event data from FITS files.
It handles the complexities of different file formats and provides a consistent
interface for accessing event data.

# Arguments
- `path::AbstractString`: Path to the FITS file

# Keyword Arguments
- `hdu::Int = 2`: HDU index to read from (typically 2 for event data)
- `T::Type = Float64`: Type to cast the time and energy columns to
- `sort::Bool = false`: Whether to sort by time if not already sorted
- `extra_columns::Vector{String} = []`: Extra columns to read from the same HDU
- `energy_alternatives::Vector{String} = ["ENERGY", "PI", "PHA"]`: Energy column alternatives to try

# Returns
`EventList{Vector{T}, FITSMetadata{FITSIO.FITSHeader}}`: EventList containing the event data

# Examples
```julia
# Basic usage
ev = readevents("events.fits")

# With custom options
ev = readevents("events.fits", hdu=3, sort=true, T=Float32)

# Reading extra columns
ev = readevents("events.fits", extra_columns=["RAWX", "RAWY", "DETX", "DETY"])

# Accessing the data
println("Number of events: \$(length(ev))")
println("Time range: \$(extrema(times(ev)))")
if has_energies(ev)
    println("Energy range: \$(extrema(energies(ev)))")
    println("Energy column: \$(ev.meta.energy_units)")
end
```

# Type Stability
This function is designed to be type-stable with proper type annotations
on return values from FITS reading operations. The return type is fully
specified to enable compiler optimizations.

# Error Handling
- Throws `AssertionError` if time and energy vectors have different sizes
- Throws `AssertionError` if times are not sorted and `sort=false`
- FITS reading errors are propagated from the FITSIO.jl library

# Implementation Notes
- Uses type-stable FITS reading with explicit type conversions
- Handles missing energy data gracefully
- Supports efficient multi-column sorting when `sort=true`
- Creates metadata with all relevant file information
- Validates data consistency before returning
- Added case_sensitive=false parameter: This tells FITSIO.jl to use the old behavior for backward compatibility
"""
function readevents(
    path::AbstractString;
    hdu::Int = 2,
    T::Type = Float64,
    sort::Bool = false,
    extra_columns::Vector{String} = String[],
    energy_alternatives::Vector{String} = ["ENERGY", "PI", "PHA"],
    kwargs...,
)::EventList{Vector{T},FITSMetadata{FITSIO.FITSHeader}}

    # Read data from FITS file with type-stable operations
    time::Vector{T},
    energy::Union{Nothing,Vector{T}},
    energy_col::Union{Nothing,String},
    header::FITSIO.FITSHeader,
    extra_data::Dict{String,Vector} = FITS(path, "r") do f

        selected_hdu = f[hdu]

        # Read header (type-stable)
        header = read_header(selected_hdu)

        # Get actual column names to find the correct TIME column
        all_cols = FITSIO.colnames(selected_hdu)
        time_col = findfirst(col -> uppercase(col) == "TIME", all_cols)
        
        if isnothing(time_col)
            error("TIME column not found in HDU $hdu")
        end
        
        actual_time_col = all_cols[time_col]
        
        # Read time column with case-insensitive option
        time = convert(Vector{T}, read(selected_hdu, actual_time_col, case_sensitive=false))

        # Read energy column using separated function
        energy_column, energy = read_energy_column(
            selected_hdu;
            T = T,
            energy_alternatives = energy_alternatives,
        )

        # Read extra columns with case-insensitive option
        extra_data = Dict{String,Vector}()
        for col_name in extra_columns
            # Find actual column name (case-insensitive)
            actual_col_idx = findfirst(col -> uppercase(col) == uppercase(col_name), all_cols)
            if !isnothing(actual_col_idx)
                actual_col_name = all_cols[actual_col_idx]
                extra_data[col_name] = read(selected_hdu, actual_col_name, case_sensitive=false)
            else
                @warn "Column '$col_name' not found in FITS file"
            end
        end

        (time, energy, energy_column, header, extra_data)
    end

    # Validate energy-time consistency
    if !isnothing(energy)
        @assert size(time) == size(energy) "Time and energy do not match sizes ($(size(time)) != $(size(energy)))"
    end

    # Handle sorting if requested
    if !issorted(time)
        if sort
            # Efficient sorting of multiple arrays
            sort_indices = sortperm(time)
            time = time[sort_indices]

            if !isnothing(energy)
                energy = energy[sort_indices]
            end

            # Sort extra columns
            for (col_name, col_data) in extra_data
                extra_data[col_name] = col_data[sort_indices]
            end
        else
            @assert false "Times are not sorted (pass `sort = true` to force sorting)"
        end
    end

    # Create metadata - just record the column name that was found
    meta = FITSMetadata(path, hdu, energy_col, extra_data, header)

    # Return type-stable EventList
    EventList(time, energy, meta)
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    summary(ev::EventList)

Provide a comprehensive summary of the EventList contents.

# Arguments
- `ev::EventList`: EventList to summarize

# Returns
String with summary information including:
- Number of events
- Time span
- Energy range (if available)
- Energy units (if available)
- Number of extra columns

# Examples
```julia
ev = readevents("events.fits")
println(summary(ev))
# Output: "EventList: 1000 events over 3600.0 time units, energies: 0.5 - 12.0 (PI), 2 extra columns"
```
"""
function Base.summary(ev::EventList)
    n_events = length(ev)
    time_span = isempty(ev.times) ? 0.0 : maximum(ev.times) - minimum(ev.times)

    summary_str = "EventList: $n_events events over $(time_span) time units"

    if has_energies(ev)
        energy_range = extrema(ev.energies)
        summary_str *= ", energies: $(energy_range[1]) - $(energy_range[2])"
        if !isnothing(ev.meta.energy_units)
            summary_str *= " ($(ev.meta.energy_units))"
        end
    end

    if !isempty(ev.meta.extra_columns)
        summary_str *= ", $(length(ev.meta.extra_columns)) extra columns"
    end

    return summary_str
end
