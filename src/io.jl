"""
    load_events_from_fits(filename::String)
Loads event list from HDU EVENTS of file fits_file, with Good Time
intervals.
"""
function load_events_from_fits(filename::String)
    FITS(filename) do hduList
        eventHDU = hduList["EVENTS"]
        cols = FITSIO.colnames(eventHDU)
        df = DataFrame(eventHDU)
        if "TIME" in cols
            time = df[!,"TIME"]
        else
            throw(KeyError("Time stamp data not provided or column names not appropriate in the file"))
        end
        ev = EventList(time=time)
        if "PI" in cols
            PI = df[!,"PI"]
            ev.PI = PI
        end
        if "ENERGY" in cols
            energy = df[!,"ENERGY"]
            ev.energy = energy
        end
        ev.gti = load_gtis(filename)
        header = FITSIO.read_header(eventHDU)
        if haskey(header,"MJDREFI") 
            ev.mjdref += header["MJDREFI"]
        end
        if haskey(header,"MJDREFF") 
            ev.mjdref += header["MJDREFF"]
        end
        if haskey(header,"INSTRUME") 
            ev.instr = header["INSTRUME"]
        end
        if haskey(header,"TIMEREF") 
            ev.timeref = header["TIMEREF"]
        end
        if haskey(header,"TIMESYS") 
            ev.timesys = header["TIMESYS"]
        end
        if haskey(header,"PLEPHEM") 
            ev.ephem = header["PLEPHEM"]
        end
        ev
    end
end

"""
    write_events_to_fits(filename::String, ev::EventList)
Writes Event and GTI HDUs to a file from provided EventList
"""
function write_events_to_fits(filename::String, ev::EventList)
    FITS(filename,"w") do hduList
        
        tkeys = String[]
        values = []

        for field in fieldnames(EventList)
            fval = getfield(ev, field)
            if !(fval isa AbstractVecOrMat)
                push!(tkeys, String(field))
                push!(values, fval)
            end
        end

        header = FITSHeader(tkeys, values,fill("",length(tkeys)))

        eventData = Dict("TIME"=>ev.time, "ENERGY"=>ev.energy, "PI"=>ev.PI);
        FITSIO.write(hduList, eventData, header = header, name = "EVENTS")

        gtiData = Dict("START"=>ev.gti[:,1], "STOP"=>ev.gti[:,2])
        FITSIO.write(hduList, gtiData, name = "GTI")
    end
end
