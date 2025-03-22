struct Meta
    headers::Vector{Dict{String,Any}}
end

struct EventList{T}
    filename::String
    times::Vector{T}
    energies::Vector{T}
    metadata::Meta
end

function readevents(path; T = Float64)
    headers = Dict{String,Any}[]
    times = T[]
    energies = T[]

    FITS(path, "r") do f
        for i = 1:length(f)  # Iterate over HDUs
            hdu = f[i]
            header_dict = Dict{String,Any}()
            for key in keys(read_header(hdu))
                header_dict[string(key)] = read_header(hdu)[key]
            end
            push!(headers, header_dict)

            # Check if the HDU is a table
            if isa(hdu, TableHDU)
                # Get column names using the correct FITSIO method
                colnames = FITSIO.colnames(hdu)

                if "TIME" in colnames
                    times = convert(Vector{T}, read(hdu, "TIME"))
                end
                if "ENERGY" in colnames
                    energies = convert(Vector{T}, read(hdu, "ENERGY"))
                end
            end
        end
    end

    metadata = Meta(headers)
    return EventList{T}(path, times, energies, metadata)
end
