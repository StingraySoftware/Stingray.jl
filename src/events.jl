"""
    FITSMetadata{H}

Metadata associated with a FITS or events file (`.fits` or `.evt`).

$(FIELDS)
"""
struct FITSMetadata{H}
    "Path to the FITS file."
    filepath::String
    "HDU index that the metadata was read from."
    hdu::Int
    "Units of energy (currently just ENERGY or PI or HPA)"
    energy_units::Union{Nothing,String}
    "Extra columns that were requested during read."
    extra_columns::Dict{String,Vector}
    "Fits headers from the selected HDU."
    headers::H
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(m::FITSMetadata))
    println(
        io,
        "FITSMetadata for $(m.filepath)[$(m.hdu)] with $(length(m.extra_columns)) extra column(s)",
    )
end

"""
    EventList{TimeType, MetaType <: FITSMetadata}

Container for an events list. Generally should not be directly constructed, but
read from file using [`readevents`](@ref).

$(FIELDS)
"""
struct EventList{TimeType<:AbstractVector,MetaType<:FITSMetadata}
    "Vector with recorded times."
    times::TimeType
    "Vector with recorded energies (else `nothing`)"
    energies::Union{Nothing,TimeType}
    meta::MetaType
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(ev::EventList))
    print(io, "EventList with $(length(ev.times)) times")
    if !isnothing(ev.energies)
        print(io, " and energies")
    end
    println(io)
end

"""
    Base.filter!(f, ev::EventList)

Filter all columns of the eventlist based on a predicate `f` applied to the
times.

# Example

```julia
# filter only positive times
filter!(t -> t > 0, ev)
```
"""
function Base.filter!(f, ev::EventList)
    # modified from the Base.filter! implementation
    j = firstindex(ev.times)

    for i in eachindex(ev.times)
        ev.times[j] = ev.times[i]

        if !isnothing(ev.energies)
            ev.energies[j] = ev.energies[i]
        end

        for (_, col) in ev.meta.extra_columns
            col[j] = col[i]
        end

        j = ifelse(f(ev.times[i])::Bool, nextind(ev.times, j), j)
    end

    if j <= lastindex(ev.times)
        resize!(ev.times, j-1)

        if !isnothing(ev.energies)
            resize!(ev.energies, j-1)
        end

        for (_, col) in ev.meta.extra_columns
            resize!(col, j-1)
        end
    end

    ev
end

"""
    colnames(file::AbstractString; hdu = 2)

Return a vector of all column names of `file`, reading a specific data unit
`hdu` if applicable.
"""
function colnames(file::AbstractString; hdu = 2)
    FITS(file) do f
        selected_hdu = f[hdu]
        FITSIO.colnames(selected_hdu)
    end
end

"""
    read_energy_column(hdu; energy_alternatives = ["ENERGY", "PI", "PHA"], T = Float64)

Attempt to read the energy column of a HDU from a list of alternative names.
Returns `Tuple{String,Vector{T}}` with the column that applied, and the column
cast to `Vector{T}`.

If no column applies, returns `(nothing, nothing)`.
"""
function read_energy_column(hdu; energy_alternatives = ["ENERGY", "PI", "PHA"], T = Float64)
    all_cols = uppercase.(FITSIO.colnames(hdu))
    for col_name in energy_alternatives
        if col_name in all_cols
            return col_name, convert.(T, read(hdu, col_name))
        end
    end
    nothing, nothing
end

"""
    readevents(path; kwargs...)

Read an [`EventList`](@ref) from a FITS file. Will attempt to read an energy
column if one exists.

Keyword arguments and defaults:
- `hdu = 2`: which HDU unit to read
- `T = Float64`: the type to cast the time and energy columns to.
- `extra_columns = []`: extra columns to read from the same HDU.

All other `kwargs` are passed to [`Stingray.read_energy_column`](@ref).
"""
function readevents(
    path::AbstractString;
    hdu = 2,
    T = Float64,
    sort = false,
    extra_columns = [],
    kwargs...,
)
    time, energy, ecol, header, add = FITS(path, "r") do f
        selected_hdu = f[hdu]
        header = read_header(selected_hdu)
        time = convert.(T, read(selected_hdu, "TIME"))
        energy_column, energy = read_energy_column(selected_hdu; T = T, kwargs...)

        add = Dict{String,Vector}(
            name => read(selected_hdu, name) for name in extra_columns
        )

        (time, energy, energy_column, header, add)
    end

    if !isnothing(energy)
        @assert size(time) == size(energy) "Time and energy do not match sizes ($(size(time)) != $(size(energy)))."
    end

    if !issorted(time)
        if sort
            # cheap and cheerful way to sort two arrays
            # there are more performant ways to do this if it becomes an issue
            I = sortperm(time)
            time = time[I]
            if !isnothing(energy)
                energy = energy[I]
            end
            for (k, col) in add
                add[k] = col[I]
            end
        else
            @assert false "Times are not sorted (pass `sort = true` to force sorting)."
        end
    end

    EventList(
        time::Vector{T},
        energy,
        FITSMetadata(path, hdu, ecol, add, header::FITSIO.FITSHeader),
    )
end
