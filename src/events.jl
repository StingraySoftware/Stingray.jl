import Base: read, write, join, sort, sort!

"""
    EventList(time::AbstractVector)
    EventList(time::T1, energy::T2=Float64[], ncounts::Int=0, mjdref::Float64=0,
              dt::Float64=0, notes::String="", gti::T3=reshape(Float64[],0,2)
              PI::Vector{Int}=Int[], high_precision::Bool=false, mission::String=""
              instr::String="", header::String="", detector_id::Vector{String}=String[]
              ephem::String="", timeref::String="", timesys::String="") where 
              {T1<:AbstractVector,T2<:AbstractVector,T3<:AbstractMatrix}
Makes an event list object. Event lists generally correspond to individual events (e.g. photons)
recorded by the detector, and their associated properties. For X-ray data where this type commonly occurs,
events are time stamps of when a photon arrived in the detector, and (optionally) the photon energy associated
with the event
"""
@with_kw mutable struct EventList{T1<:AbstractVector,T2<:AbstractVector,T3<:AbstractMatrix}
    time::T1
    energy::T2=Float64[]
    ncounts::Int=0
    mjdref::Float64=0
    dt::Float64=0
    notes::String=""
    gti::T3=reshape(Float64[],0,2)
    PI::Vector{Int}=Int[]
    high_precision::Bool=false
    mission::String=""
    instr::String=""
    header::String=""
    detector_id::Vector{String}=String[]
    ephem::String=""
    timeref::String=""
    timesys::String=""
end

"""
    read(::Type{EventList},filename::String, format::String)
Read the EventList struct from a file. Only FITS file is supported currently.
"""
function read(::Type{EventList},filename::String, format::String)
    if format=="fits"
        return load_events_from_fits(filename)
    else
        throw(ArgumentError("File format $(format) not supported still."))
    end
end

"""
    write(ev::EventList, filename::String, format::String)
Write EventList to a file following its format. Only FITS file is supported currently.
"""
function write(ev::EventList, filename::String, format::String)
    if format=="fits"
        write_events_to_fits(filename, ev)
    else
        throw(ArgumentError("File format $(format) not supported still."))
    end
end

"""
    from_lc(lc::LightCurve)
Create an `EventList` object from a `Lightcurve` object. 
"""
function from_lc(lc::LightCurve)
    times = [lc.time[i] for i in 1:length(lc.time) for _ in 1:lc.counts[i]]
    return EventList(time=times,gti=lc.gti)
end

"""
    to_lc(ev::EventList, dt)
    to_lc(ev::EventList, dt::Real, tstart::Real, tseg::Real)
Convert `EventList` object to a `Lightcurve` object.
"""
function to_lc(ev::EventList, dt)
    if isempty(ev.gti)
        tstart=0
        tseg=0
    else
        tstart = ev.gti[begin][begin]
        tseg = ev.gti[end][end]-tstart
    end
    return to_lc(ev, dt, tstart, tseg)
end

to_lc(ev::EventList, dt::Real,
      tstart::Real, tseg::Real) = make_lightcurves(ev.time, dt, ev.gti;
                                                   tstart=tstart, tseg=tseg,
                                                   mjdref=ev.mjdref)

"""
    join(ev1::EventList, ev2::EventList)
Join two `EventList` objects into one.
If both are empty, an empty `EventList` is returned.
GTIs are crossed if the event lists are over a common time interval,
and appended otherwise.
"""
function join(ev1::EventList, ev2::EventList)
    new_ev = EventList(time=vcat(ev1.time, ev2.time))
    if ev1.dt!=ev2.dt
        new_ev.dt = max(ev1.dt,ev2.dt);
    end

    # Tolerance for MJDREF:1 microsecond
    if !isapprox(ev1.mjdref, ev2.mjdref, atol=1e-6 / 86400)
        ev2.mjdref = ev1.mjdref
    end

    order = sortperm(new_ev.time)
    new_ev.time = new_ev.time[order]

    if !isempty(ev1.PI) || !isempty(ev2.PI)
        if isempty(ev1.PI)
            ev1.PI = zero(ev1.time)
        end
        if isempty(ev2.PI)
            ev2.PI = zero(ev2.time)
        end
        new_ev.PI = vcat(ev1.PI, ev2.PI)
        new_ev.PI = new_ev.PI[order]
    else
        new_ev.PI = Float64[]
    end

    if !isempty(ev1.energy) || !isempty(ev2.energy)
        if isempty(ev1.energy)
            ev1.energy = zero(ev1.time)
        end
        if isempty(ev2.energy)
            ev2.energy = zero(ev2.time)
        end
        new_ev.energy = vcat(ev1.energy, ev2.energy)
        new_ev.energy = new_ev.energy[order]
    else
        new_ev.energy = Float64[]
    end

    if !isempty(ev1.gti) || !isempty(ev2.gti)
        if isempty(ev1.gti)
            ev1.gti = [ev1.time[begin] - ev1.dt/2 ev1.time[end] + ev1.dt/2]
        end
        if isempty(ev2.gti)
            ev2.gti = [ev2.time[begin] - ev2.dt/2 ev2.time[end] + ev2.dt/2]
        end
        new_ev.gti = operations_on_gtis([ev1.gti,ev2.gti], intersect)
        if isempty(new_ev.gti)
            @warn """
            GTIs in these two event lists do not overlap at all.
            Merging instead of returning an overlap.
            """
            new_ev.gti = operations_on_gtis([ev1.gti,ev2.gti], union)
        end
    else
        new_ev.gti = reshape(Float64[],0,2)
    end

    for attr in (:mission, :instr)
        if getfield(ev1, attr) != getfield(ev2, attr)
            setfield!(new_ev, attr, string(getfield(ev1, attr),',',getfield(ev2, attr)))
        else
            setfield!(new_ev, attr, getfield(ev1, attr))
        end
    end

    new_ev.mjdref = ev1.mjdref

    return new_ev
end
        
"""
    sort(ev::EventList)
Sort the event list in time and return a new event list.
"""
function sort(ev::EventList)
    order = sortperm(ev.time)
    new_ev = deepcopy(ev)
    apply_mask(new_ev, order)
    return new_ev
end

"""
    sort!(ev::EventList)
Sort (in place) the event list in time.
"""
function sort!(ev::EventList)
    order = sortperm(ev.time)
    apply_mask(ev, order)
end

"""
    apply_mask(ev::EventList, mask::AbstractVector{T};) where {T<:Union{Bool,Int}}
Apply a mask to all array attributes of the event list
"""
function apply_mask(ev::EventList, mask::AbstractVector{T};) where {T<:Union{Bool,Int}}
    for field in fieldnames(EventList)
        fval = getfield(ev, field)
        if fval isa AbstractVector && !isempty(fval) 
            setfield!(ev, field, deepcopy(fval)[mask])
        end
    end
end
