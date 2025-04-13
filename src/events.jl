"""
Abstract type representing a power spectrum, which characterizes the distribution 
of power across different frequencies in a signal.

Subtypes include:
- PowerSpectrum{T}: Represents a power spectrum for a single signal segment
- AveragedPowerspectrum{T}: Represents a power spectrum averaged over multiple segments
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
- `energies::Vector{T}`: Vector of event energies.
- `metadata::DictMetadata`: Metadata information extracted from the FITS file headers.
"""
struct EventList{T} <: AbstractEventList{T}
    filename::String
    times::Vector{T}
    energies::Vector{T}
    metadata::DictMetadata
end

times(ev::EventList) = ev.times
energies(ev::EventList) = ev.energies

"""
    readevents(path; T = Float64)

Read event data from a FITS file into an EventList structure.

## Arguments
- `path::String`: Path to the FITS file
- `T::Type=Float64`: Numeric type for the data

## Returns
- [`EventList`](@ref) containing the extracted data

## Notes

The function extracts `TIME` and `ENERGY` columns from any TableHDU in the FITS
file. All headers from each HDU are collected into the metadata field. It will
use the first occurrence of complete event data (both TIME and ENERGY columns)
found in the file.
"""
function readevents(path; T = Float64)
    headers = Dict{String,Any}[]
    times = T[]
    energies = T[]

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
                end
                if isempty(energies) && ("ENERGY" in colnames)
                    energies = convert(Vector{T}, read(hdu, "ENERGY"))
                end

                # If we found both time and energy data, we can return
                if !isempty(times) && !isempty(energies)
                    @info "Found complete event data in extension $(i) of $(path)"
                    metadata = DictMetadata(headers)
                    return EventList{T}(path, times, energies, metadata)
                end
            end
        end
    end

    if isempty(times)
        @warn "No TIME data found in FITS file $(path). Time series analysis will not be possible."
    end
    if isempty(energies)
        @warn "No ENERGY data found in FITS file $(path). Energy spectrum analysis will not be possible."
    end

    metadata = DictMetadata(headers)
    return EventList{T}(path, times, energies, metadata)
end

Base.length(ev::AbstractEventList) = length(times(ev))
Base.size(ev::AbstractEventList) = (length(ev),)
Base.getindex(ev::EventList, i) = (ev.times[i], ev.energies[i])

function Base.show(io::IO, ev::EventList{T}) where T
    print(io, "EventList{$T}(n=$(length(ev)), file=$(ev.filename))")
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
