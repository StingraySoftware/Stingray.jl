using FITSIO

"""
Abstract type for all event list implementations
"""
abstract type AbstractEventList{T} end

"""
    DictMetadata

A structure containing metadata from FITS file headers.

## Fields

- `headers::Vector{Dict{String,Any}}`: A vector of dictionaries containing header information from each HDU.
"""
struct DictMetadata
    headers::Vector{Dict{String,Any}}
end

"""
    EventList{T} <: AbstractEventList{T}

A structure containing event data from a FITS file.

## Fields

- `filename::String`: Path to the source FITS file.
- `times::Vector{T}`: Vector of event times.
- `energies::Union{Vector{T}, Nothing}`: Vector of event energies (or nothing if not available).
- `extra_columns::Dict{String, Vector}`: Dictionary of additional column data.
- `metadata::DictMetadata`: Metadata information extracted from the FITS file headers.
"""
struct EventList{T} <: AbstractEventList{T}
    filename::String
    times::Vector{T}
    energies::Union{Vector{T}, Nothing}
    extra_columns::Dict{String, Vector}
    metadata::DictMetadata
end

# Simplified constructor that defaults energies to nothing and extra_columns to empty dict
function EventList{T}(filename, times, metadata) where T
    EventList{T}(filename, times, nothing, Dict{String, Vector}(), metadata)
end

# Constructor that accepts energies but defaults extra_columns to empty dict
function EventList{T}(filename, times, energies, metadata) where T
    EventList{T}(filename, times, energies, Dict{String, Vector}(), metadata)
end

times(ev::EventList) = ev.times
energies(ev::EventList) = ev.energies

"""
    readevents(path; T = Float64, energy_alternatives=["ENERGY", "PI", "PHA"])

Read event data from a FITS file into an EventList structure.

## Arguments
- `path::String`: Path to the FITS file
- `T::Type=Float64`: Numeric type for the data
- `energy_alternatives::Vector{String}=["ENERGY", "PI", "PHA"]`: Column names to try for energy data

## Returns
- [`EventList`](@ref) containing the extracted data

## Notes

The function extracts `TIME` and energy columns from TableHDUs in the FITS
file. All headers from each HDU are collected into the metadata field. When it finds
an HDU containing TIME column, it also looks for energy data and collects additional
columns from that same HDU, since all event data is typically stored together.
"""
function readevents(path; T = Float64, energy_alternatives=["ENERGY", "PI", "PHA"])
    headers = Dict{String,Any}[]
    times = T[]
    energies = T[]
    extra_columns = Dict{String, Vector}()
    
    FITS(path, "r") do f
        for i = 1:length(f)  # Iterate over HDUs
            hdu = f[i]
            
            # Always collect headers from all extensions
            header_dict = Dict{String,Any}()
            for key in keys(read_header(hdu))
                header_dict[string(key)] = read_header(hdu)[key]
            end
            push!(headers, header_dict)
            
            # Check if the HDU is a table
            if isa(hdu, TableHDU)
                colnames = FITSIO.colnames(hdu)
                
                # Read TIME and ENERGY data if columns exist and vectors are empty
                if isempty(times) && ("TIME" in colnames)
                    times = convert(Vector{T}, read(hdu, "TIME"))
                    @info "Found TIME column in extension $(i) of $(path)"
                    
                    # Once we find the TIME column, process all other columns in this HDU
                    # as this is where all event data will be
                    
                    # Try ENERGY first
                    if "ENERGY" in colnames && isempty(energies)
                        energies = convert(Vector{T}, read(hdu, "ENERGY"))
                        @info "Found ENERGY column in the same extension"
                    else
                        # Try alternative energy columns if ENERGY is not available
                        for energy_col in energy_alternatives[2:end]  # Skip ENERGY as we already checked
                            if energy_col in colnames && isempty(energies)
                                energies = convert(Vector{T}, read(hdu, energy_col))
                                @info "Using '$energy_col' column for energy information"
                                break
                            end
                        end
                    end
                    
                    # Collect all columns from this HDU for extra_columns
                    for col in colnames
                        # Add every column to extra_columns for consistent access
                        try
                            extra_columns[col] = read(hdu, col)
                            @debug "Added column '$col' to extra_columns"
                        catch e
                            @warn "Failed to read column '$col': $e"
                        end
                    end
                    
                    # We've found and processed the event data HDU, stop searching
                    break
                end
            end
        end
    end
    
    if isempty(times)
        @warn "No TIME data found in FITS file $(path). Time series analysis will not be possible."
    end
    if isempty(energies)
        @warn "No ENERGY data found in FITS file $(path). Energy spectrum analysis will not be possible."
        energies = nothing
    end
    
    metadata = DictMetadata(headers)
    return EventList{T}(path, times, energies, extra_columns, metadata)
end


Base.length(ev::AbstractEventList) = length(times(ev))
Base.size(ev::AbstractEventList) = (length(ev),)

function Base.getindex(ev::EventList, i)
    if isnothing(ev.energies)
        return (ev.times[i], nothing)
    else
        return (ev.times[i], ev.energies[i])
    end
end

function Base.show(io::IO, ev::EventList{T}) where T
    energy_status = isnothing(ev.energies) ? "no energy data" : "with energy data"
    extra_cols = length(keys(ev.extra_columns))
    print(io, "EventList{$T}(n=$(length(ev)), $energy_status, $extra_cols extra columns, file=$(ev.filename))")
end

"""
    validate(events::AbstractEventList)

Validate the event list structure.

## Returns
- `true` if valid, throws ArgumentError otherwise
"""
function validate(events::AbstractEventList)
    evt_times = times(events)
    if !issorted(evt_times)
        throw(ArgumentError("Event times must be sorted in ascending order"))
    end
    if length(evt_times) == 0
        throw(ArgumentError("Event list is empty"))
    end
    return true
end


"""
    get_column(events::EventList, column_name::String)

Get a specific column from the event list.

## Arguments
- `events::EventList`: Event list object
- `column_name::String`: Name of the column to retrieve

## Returns
- The column data if available, nothing otherwise
"""
function get_column(events::EventList, column_name::String)
    if column_name == "TIME"
        return events.times
    elseif column_name == "ENERGY" && !isnothing(events.energies)
        return events.energies
    elseif haskey(events.extra_columns, column_name)
        return events.extra_columns[column_name]
    else
        return nothing
    end
end