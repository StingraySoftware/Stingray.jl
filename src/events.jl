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
    EventList{T}

A structure containing event data from a FITS file.

## Fields

- `filename::String`: Path to the source FITS file.
- `times::Vector{T}`: Vector of event times.
- `energies::Vector{T}`: Vector of event energies.
- `metadata::DictMetadata`: Metadata information extracted from the FITS file headers.
"""
struct EventList{T}
    filename::String
    times::Vector{T}
    energies::Vector{T}
    metadata::DictMetadata
end

"""
    readevents(path; T = Float64)

Read event data from a FITS file into an EventList structure. The `path` is a
string that points to the location of the FITS file. `T` is used to specify
which numeric type to convert the data to.

Returns an [`EventList`](@ref) containing the extracted data.

## Notes

The function extracts `TIME` and `ENERGY` columns from any TableHDU in the FITS
file. All headers from each HDU are collected into the metadata field.
"""
function readevents(path; T = Float64)
    headers = Dict{String,Any}[]  # Store headers from all HDUs
    times = T[]  # Store times
    energies = T[]  # Store energies
    event_data_read = false  # Flag to check if event data has already been read

    FITS(path, "r") do f
        for i = 1:length(f)  # Iterate over all HDUs
            hdu = f[i]
            header_dict = Dict{String,Any}()
            for key in keys(read_header(hdu))
                header_dict[string(key)] = read_header(hdu)[key]
            end
            push!(headers, header_dict)  # Collect header information for each HDU

            # Check if the HDU is a table
            if isa(hdu, TableHDU)
                colnames = FITSIO.colnames(hdu)

                # If both TIME and ENERGY columns are found and we haven't read the event data yet
                if "TIME" in colnames && "ENERGY" in colnames && !event_data_read
                    times = convert(Vector{T}, read(hdu, "TIME"))
                    energies = convert(Vector{T}, read(hdu, "ENERGY"))
                    event_data_read = true  # Mark that we've read the event data
                elseif "TIME" in colnames && !event_data_read
                    times = convert(Vector{T}, read(hdu, "TIME"))
                elseif "ENERGY" in colnames && !event_data_read
                    energies = convert(Vector{T}, read(hdu, "ENERGY"))
                end
            end
        end
    end

    # Check if no data found for times or energies
    if isempty(times)
        @warn "No TIME data found in FITS file $(path). Time series analysis will not be possible."
    end

    if isempty(energies)
        @warn "No ENERGY data found in FITS file $(path). Energy spectrum analysis will not be possible."
    end

    metadata = DictMetadata(headers)  # Metadata for all headers
    return EventList{T}(path, times, energies, metadata)
end

