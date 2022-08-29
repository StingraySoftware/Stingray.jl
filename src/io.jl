function load_events_from_fits(filename::String)
    a = FITS(filename) do hduList
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
        if "MJDREFI" in header.keys
            ev.mjdref += header["MJDREFI"]
        end
        if "MJDREFF" in header.keys
            ev.mjdref += header["MJDREFF"]
        end
        if "INSTRUME" in header.keys
            ev.instr = header["INSTRUME"]
        end
        if "TIMEREF" in header.keys
            ev.timeref = header["TIMEREF"]
        end
        if "TIMESYS" in header.keys
            ev.timesys = header["TIMESYS"]
        end
        if "PLEPHEM" in header.keys
            ev.ephem = header["PLEPHEM"]
        end
        ev
    end
    return a
end

# TODO
function write_events_to_fits(filename::String, ev::EventList)
    # FITS(filename,"w") do event
    #      FITSIO.write(event, ev)
    # end
end