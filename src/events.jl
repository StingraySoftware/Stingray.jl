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

    # Inner constructor with validation
    function EventList{T}(filename::String, times::Vector{T}, energies::Union{Vector{T}, Nothing}, 
                         extra_columns::Dict{String, Vector}, metadata::DictMetadata) where T
        # Validate event times
        if isempty(times)
            throw(ArgumentError("Event list cannot be empty"))
        end
        
        if !issorted(times)
            throw(ArgumentError("Event times must be sorted in ascending order"))
        end
        
        # Validate energy vector length if present
        if !isnothing(energies) && length(energies) != length(times)
            throw(ArgumentError("Energy vector length ($(length(energies))) must match times vector length ($(length(times)))"))
        end
        
        # Validate extra columns have consistent lengths
        for (col_name, col_data) in extra_columns
            if length(col_data) != length(times)
                throw(ArgumentError("Column '$col_name' length ($(length(col_data))) must match times vector length ($(length(times)))"))
            end
        end
        
        new{T}(filename, times, energies, extra_columns, metadata)
    end
end

# Simplified constructors that use the validated inner constructor
function EventList{T}(filename, times, metadata) where T
    EventList{T}(filename, times, nothing, Dict{String, Vector}(), metadata)
end

function EventList{T}(filename, times, energies, metadata) where T
    EventList{T}(filename, times, energies, Dict{String, Vector}(), metadata)
end

# Accessor functions
times(ev::EventList) = ev.times
energies(ev::EventList) = ev.energies

"""
    readevents(path; T = Float64, energy_alternatives=["ENERGY", "PI", "PHA"])

Read event data from a FITS file into an EventList structure with enhanced performance.

## Arguments
- `path::String`: Path to the FITS file
- `T::Type=Float64`: Numeric type for the data
- `energy_alternatives::Vector{String}=["ENERGY", "PI", "PHA"]`: Column names to try for energy data

## Returns
- [`EventList`](@ref) containing the extracted data
"""
function readevents(path::String;
                   mission::Union{String,Nothing}=nothing,
                   instrument::Union{String,Nothing}=nothing,
                   epoch::Union{Float64,Nothing}=nothing,
                   T::Type=Float64,
                   energy_alternatives::Vector{String}=["ENERGY", "PI", "PHA"],
                   sector_column::Union{String,Nothing}=nothing,
                   event_hdu::Int=2)  #X-ray event files have events in HDU 2
    
    # Get mission support if specified
    mission_support = if !isnothing(mission)
        ms = get_mission_support(mission, instrument, epoch)
        # Use mission-specific energy alternatives if available
        energy_alternatives = ms.energy_alternatives
        ms
    else
        nothing
    end
    
    # Initialize containers
    headers = Dict{String,Any}[]
    times = T[]
    energies = T[]
    extra_columns = Dict{String, Vector}()
    
    FITS(path, "r") do f
        # Collect headers from all HDUs
        for i = 1:length(f)
            hdu = f[i]
            header_dict = Dict{String,Any}()
            try
                for key in keys(read_header(hdu))
                    header_dict[string(key)] = read_header(hdu)[key]
                end
            catch e
                @debug "Could not read header from HDU $i: $e"
            end
            
            # Apply mission-specific patches to header information
            if !isnothing(mission)
                header_dict = patch_mission_info(header_dict, mission)
            end
            push!(headers, header_dict)
        end
        
        # Try to read event data from the specified HDU (default: HDU 2)
        try
            hdu = f[event_hdu]
            if !isa(hdu, TableHDU)
                throw(ArgumentError("HDU $event_hdu is not a table HDU"))
            end
            
            colnames = FITSIO.colnames(hdu)
            @info "Reading events from HDU $event_hdu with columns: $(join(colnames, ", "))"
            
            # Read TIME column (case-insensitive search)
            time_col = nothing
            for col in colnames
                if uppercase(col) == "TIME"
                    time_col = col
                    break
                end
            end
            
            if isnothing(time_col)
                throw(ArgumentError("No TIME column found in HDU $event_hdu"))
            end
            
            # Read time data
            raw_times = read(hdu, time_col)
            times = convert(Vector{T}, raw_times)
            @info "Successfully read $(length(times)) events"
            
            # Try to read energy data
            energy_col = nothing
            for ecol in energy_alternatives
                for col in colnames
                    if uppercase(col) == uppercase(ecol)
                        energy_col = col
                        @info "Using '$col' column for energy data"
                        break
                    end
                end
                if !isnothing(energy_col)
                    break
                end
            end
            
            if !isnothing(energy_col)
                try
                    raw_energy = read(hdu, energy_col)
                    energies = if !isnothing(mission_support)
                        @info "Applying mission calibration for $mission"
                        convert(Vector{T}, apply_calibration(mission_support, raw_energy))
                    else
                        convert(Vector{T}, raw_energy)
                    end
                    @info "Energy data: $(length(energies)) values, range: $(extrema(energies))"
                catch e
                    @warn "Failed to read energy column '$energy_col': $e"
                    energies = T[]
                end
            else
                @info "No energy column found in available alternatives: $(join(energy_alternatives, ", "))"
            end
            
            # Read additional columns if specified
            if !isnothing(sector_column)
                sector_col_found = nothing
                for col in colnames
                    if uppercase(col) == uppercase(sector_column)
                        sector_col_found = col
                        break
                    end
                end
                
                if !isnothing(sector_col_found)
                    try
                        extra_columns["SECTOR"] = read(hdu, sector_col_found)
                        @info "Read sector/detector data from '$sector_col_found'"
                    catch e
                        @warn "Failed to read sector column '$sector_col_found': $e"
                    end
                end
            end
            
        catch e
            # If default HDU fails, fall back to searching all HDUs
            @warn "Failed to read from HDU $event_hdu: $e. Searching all HDUs..."
            
            event_found = false
            for i = 1:length(f)
                hdu = f[i]
                if isa(hdu, TableHDU)
                    try
                        colnames = FITSIO.colnames(hdu)
                        # Look for TIME column
                        if any(uppercase(col) == "TIME" for col in colnames)
                            @info "Found events in HDU $i"
                            raw_times = read(hdu, "TIME")
                            times = convert(Vector{T}, raw_times)
                            
                            # Try to read energy
                            for ecol in energy_alternatives
                                for col in colnames
                                    if uppercase(col) == uppercase(ecol)
                                        try
                                            raw_energy = read(hdu, col)
                                            energies = convert(Vector{T}, raw_energy)
                                            break
                                        catch
                                            continue
                                        end
                                    end
                                end
                                if !isempty(energies)
                                    break
                                end
                            end
                            
                            event_found = true
                            break
                        end
                    catch
                        continue
                    end
                end
            end
            
            if !event_found
                throw(ArgumentError("No TIME column found in any HDU of FITS file $(basename(path))"))
            end
        end
    end
    
    if isempty(times)
        throw(ArgumentError("No event data found in FITS file $(basename(path))"))
    end
    
    @info "Successfully loaded $(length(times)) events from $(basename(path))"
    
    # Create metadata and return EventList
    metadata = DictMetadata(headers)
    return EventList{T}(path, 
                       times, 
                       isempty(energies) ? nothing : energies,
                       extra_columns, 
                       metadata)
end

# Basic interface methods
Base.length(ev::AbstractEventList) = length(times(ev))
Base.size(ev::AbstractEventList) = (length(ev),)

function Base.show(io::IO, ev::EventList{T}) where T
    energy_status = isnothing(ev.energies) ? "no energy data" : "with energy data"
    extra_cols = length(keys(ev.extra_columns))
    print(io, "EventList{$T}(n=$(length(ev)), $energy_status, $extra_cols extra columns, file=$(ev.filename))")
end