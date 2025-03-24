"""
    DictMetadata

A structure containing metadata from FITS file headers.

Fields
------
- `headers::Vector{Dict{String,Any}}`: A vector of dictionaries containing header information from each HDU.
"""
struct DictMetadata
    headers::Vector{Dict{String,Any}}
end

"""
    EventList{T}

A structure containing event data from a FITS file.

Fields
------
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
            
            # Check if the HDU is a table and we haven't found events yet
            if isa(hdu, TableHDU)
                colnames = FITSIO.colnames(hdu)
                has_time = "TIME" in colnames
                has_energy = "ENERGY" in colnames

                if has_time
                    times = convert(Vector{T}, read(hdu, "TIME"))

                    if has_energy
                        energies = convert(Vector{T}, read(hdu, "ENERGY"))
                    end
                    # Return immediately after finding and reading event data
                    @info "Found event data in extension $(i) of $(path)"
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
