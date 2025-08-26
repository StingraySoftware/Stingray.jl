function compare_tables(table1, table2; rtol=0.001, discard=[])
    s_discard = Symbol.(discard)
    test_result = true

    for key in propertynames(table1)
        if key in s_discard
            continue
        end
        oe, oc = getproperty(table1, key), getproperty(table2, key)
        if oe isa Integer || oe isa String
            if !(oe == oc) 
                test_result = false
                break
            end
        elseif isnothing(oe)
            if !(isnothing(oc)) 
                test_result = false
                break
            end
        else
            if !(≈(oe, oc, rtol=rtol)) 
                test_result = false
                break
            end
        end
    end

    data_fields_table1 = filter(field -> !(typeof(field) == Symbol && startswith(string(field), "meta")), names(table1))
    data_fields_table2 = filter(field -> !(typeof(field) == Symbol && startswith(string(field), "meta")), names(table2))

    for field in names(table1)
        if field in discard
            continue
        end
        oe, oc = table1[!, field], table2[!, field]

        if !(≈(oe, oc, rtol=rtol)) 
            test_result = false
            break
        end
    end
    @test test_result
end

const FOURIER_DT = 1
const FOURIER_LEN = 100
const FOURIER_CTRATE = 10000
const FOURIER_SEGMENT_SIZE = 20.0

function setup_fourier_test_data(; dt=FOURIER_DT, len=FOURIER_LEN, ctrate=FOURIER_CTRATE, two_datasets=false)
    N = len ÷ dt
    dt = len / N
    times = sort(rand(Uniform(0, len), len * ctrate))
    gti = Float64[0 len]
    bins = LinRange(0, len, N + 1)
    counts = fit(Histogram, times, bins).weights
    errs = fill!(similar(counts), 1) * sqrt(ctrate)
    bin_times = (@view(bins[1:end-1]) + @view(bins[2:end])) / 2
    
    if two_datasets
        times2 = sort(rand(Uniform(0, len), len * ctrate))
        counts2 = fit(Histogram, times2, bins).weights
        errs2 = fill!(similar(counts2), 1) * sqrt(ctrate)
        return (times=times, times2=times2, gti=gti, bins=bins, counts=counts, counts2=counts2, 
                errs=errs, errs2=errs2, bin_times=bin_times, dt=dt, N=N)
    else
        return (times=times, gti=gti, bins=bins, counts=counts, errs=errs, 
                bin_times=bin_times, dt=dt, N=N)
    end
end

function test_pds_common_mean(norm::String)
    data = setup_fourier_test_data()
    
    out_comm = avg_pds_from_events(
        data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=true, silent=true, fluxes=nothing
    ).power
    out = avg_pds_from_events(
        data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=false, silent=true, fluxes=nothing
    ).power
    @test std(out_comm) ≈ std(out) rtol=0.1
end

function test_cs_common_mean(norm::String, return_auxil::Bool, fullspec::Bool)
    data = setup_fourier_test_data(two_datasets=true)
    
    out_comm = avg_cs_from_events(
        data.times, data.times2, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=true, silent=true,
        fullspec=fullspec, return_auxil=return_auxil
    ).power
    out = avg_cs_from_events(
        data.times, data.times2, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=false, silent=true,
        fullspec=fullspec, return_auxil=return_auxil
    ).power
    @test std(out_comm) ≈ std(out) rtol=0.1
end

function test_events_vs_counts_pds(use_common_mean::Bool, norm::String, include_errors::Bool=false)
    data = setup_fourier_test_data()
    
    out_ev = avg_pds_from_events(
        data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=use_common_mean, silent=true, fluxes=nothing
    )
    
    if include_errors
        out_ct = avg_pds_from_events(
            data.bin_times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
            norm=norm, use_common_mean=use_common_mean, silent=true,
            fluxes=data.counts, errors=data.errs
        )
        rtol = use_common_mean ? 0.01 : 0.1
        compare_tables(out_ev, out_ct, rtol=rtol, discard=["variance"])
    else
        out_ct = avg_pds_from_events(
            data.bin_times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
            norm=norm, use_common_mean=use_common_mean, silent=true, fluxes=data.counts
        )
        compare_tables(out_ev, out_ct)
    end
end

function test_events_vs_counts_cs(use_common_mean::Bool, norm::String, include_errors::Bool=false)
    data = setup_fourier_test_data(two_datasets=true)
    
    out_ev = avg_cs_from_events(
        data.times, data.times2, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
        norm=norm, use_common_mean=use_common_mean, silent=true
    )
    
    if include_errors
        out_ct = avg_cs_from_events(
            data.bin_times, data.bin_times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
            norm=norm, use_common_mean=use_common_mean, silent=true,
            fluxes1=data.counts, fluxes2=data.counts2, 
            errors1=data.errs, errors2=data.errs2
        )
        discard = [m for m in propertynames(out_ev) if m == :variance]
        rtol = use_common_mean ? 0.01 : 0.1
        compare_tables(out_ev, out_ct, rtol=rtol, discard=discard)
    else
        out_ct = avg_cs_from_events(
            data.bin_times, data.bin_times, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
            norm=norm, use_common_mean=use_common_mean, silent=true,
            fluxes1=data.counts, fluxes2=data.counts2
        )
        rtol = use_common_mean ? 0.01 : 0.1
        compare_tables(out_ev, out_ct, rtol=rtol)
    end
end

let
    data = setup_fourier_test_data()
    @test get_average_ctrate(data.times, data.gti, FOURIER_SEGMENT_SIZE) == FOURIER_CTRATE
    @test get_average_ctrate(data.bin_times, data.gti, FOURIER_SEGMENT_SIZE; counts=data.counts) == FOURIER_CTRATE
end

let
    data = setup_fourier_test_data()
    N = round(Int, FOURIER_SEGMENT_SIZE / data.dt)
    fts_evts = collect(get_flux_iterable_from_segments(data.times, data.gti, FOURIER_SEGMENT_SIZE, n_bin=N))
    fts_cts = collect(get_flux_iterable_from_segments(data.bin_times, data.gti, FOURIER_SEGMENT_SIZE, fluxes=data.counts))
    @test fts_evts == fts_cts
end

let
    common_ref = true
    @test_logs (:warn, r"n_ave is below 30.") error_on_averaged_cross_spectrum([4 + 1.0im], [2], [4], 29, 2, 2; common_ref=common_ref)
end

let
    common_ref = false
    @test_logs (:warn, r"n_ave is below 30.") error_on_averaged_cross_spectrum([4 + 1.0im], [2], [4], 29, 2, 2; common_ref=common_ref)
end

let
    gti = Float64[0 FOURIER_LEN]
    _times = rand(Uniform(0, 1000), 1)
    out_ev = avg_pds_from_events(_times, gti, FOURIER_SEGMENT_SIZE, FOURIER_DT, silent=true)
    @test isnothing(out_ev)
end

let
    gti = Float64[0 FOURIER_LEN]
    return_auxil = true
    _times1 = rand(Uniform(0, 1000), 1)
    _times2 = rand(Uniform(0, 1000), 1)
    out_ev = avg_cs_from_events(_times1, _times2, gti, FOURIER_SEGMENT_SIZE, FOURIER_DT, silent=true, return_auxil=return_auxil)
    @test isnothing(out_ev) 
end

let
    gti = Float64[0 FOURIER_LEN]
    return_auxil = false
    _times1 = rand(Uniform(0, 1000), 1)
    _times2 = rand(Uniform(0, 1000), 1)
    out_ev = avg_cs_from_events(_times1, _times2, gti, FOURIER_SEGMENT_SIZE, FOURIER_DT, silent=true, return_auxil=return_auxil)
    @test isnothing(out_ev) 
end

let; test_pds_common_mean("frac"); end
let; test_pds_common_mean("abs"); end
let; test_pds_common_mean("none"); end
let; test_pds_common_mean("leahy"); end

let; test_cs_common_mean("frac", true, true); end
let; test_cs_common_mean("frac", true, false); end
let; test_cs_common_mean("frac", false, true); end
let; test_cs_common_mean("frac", false, false); end

let; test_events_vs_counts_pds(true, "frac"); end
let; test_events_vs_counts_pds(false, "frac"); end
let; test_events_vs_counts_pds(true, "frac", true); end

let; test_events_vs_counts_cs(true, "frac"); end
let; test_events_vs_counts_cs(true, "frac", true); end

# Test LightCurve with vector dt (should use first element)
let
    data = setup_fourier_test_data()
    eventlist = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    
    result_eventlist = avg_pds_from_eventlist(eventlist, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    result_times = avg_pds_from_events(data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    
    @test result_eventlist.power ≈ result_times.power
    @test result_eventlist.freq ≈ result_times.freq
end

# Test EventList PDS function equivalence with different norms
let
    data = setup_fourier_test_data()
    eventlist = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    
    for norm in ["frac", "abs", "none", "leahy"]
        result_eventlist = avg_pds_from_eventlist(eventlist, FOURIER_SEGMENT_SIZE, data.dt, norm=norm, silent=true)
        result_times = avg_pds_from_events(data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt, norm=norm, silent=true)
        
        @test result_eventlist.power ≈ result_times.power
        @test result_eventlist.freq ≈ result_times.freq
    end
end

# Test EventList PDS function with use_common_mean parameter
let
    data = setup_fourier_test_data()
    eventlist = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    
    for use_common_mean in [true, false]
        result_eventlist = avg_pds_from_eventlist(eventlist, FOURIER_SEGMENT_SIZE, data.dt, use_common_mean=use_common_mean, silent=true)
        result_times = avg_pds_from_events(data.times, data.gti, FOURIER_SEGMENT_SIZE, data.dt, use_common_mean=use_common_mean, silent=true)
        
        @test result_eventlist.power ≈ result_times.power
        @test result_eventlist.freq ≈ result_times.freq
    end
end

# Test EventList cross spectrum function equivalence
let
    data = setup_fourier_test_data(two_datasets=true)
    eventlist1 = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    eventlist2 = EventList(data.times2, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    
    result_eventlists = avg_cs_from_eventlists(eventlist1, eventlist2, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    result_times = avg_cs_from_events(data.times, data.times2, data.gti, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    
    @test result_eventlists.power ≈ result_times.power
    @test result_eventlists.freq ≈ result_times.freq
end

# Test EventList cross spectrum function with different parameters
let
    data = setup_fourier_test_data(two_datasets=true)
    eventlist1 = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    eventlist2 = EventList(data.times2, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), data.gti, nothing
    ))
    
    for norm in ["frac", "abs", "none", "leahy"]
        for use_common_mean in [true, false]
            for fullspec in [true, false]
                for return_auxil in [true, false]
                    result_eventlists = avg_cs_from_eventlists(
                        eventlist1, eventlist2, FOURIER_SEGMENT_SIZE, data.dt,
                        norm=norm, use_common_mean=use_common_mean,
                        fullspec=fullspec, return_auxil=return_auxil, silent=true
                    )
                    result_times = avg_cs_from_events(
                        data.times, data.times2, data.gti, FOURIER_SEGMENT_SIZE, data.dt,
                        norm=norm, use_common_mean=use_common_mean,
                        fullspec=fullspec, return_auxil=return_auxil, silent=true
                    )
                    
                    @test result_eventlists.power ≈ result_times.power
                    @test result_eventlists.freq ≈ result_times.freq
                    
                    if return_auxil
                        @test result_eventlists.pds1 ≈ result_times.pds1
                        @test result_eventlists.pds2 ≈ result_times.pds2
                    end
                end
            end
        end
    end
end

# Test LightCurve PDS function equivalence
let
    dt = 1.0
    N = 100
    times = collect(0:dt:(N-1)*dt)
    counts = rand(Poisson(100), N)
    
    metadata = LightCurveMetadata(
        "TEST", "INST", "TEST_OBJ", 58000.0, 
        (times[1], times[end]), dt,
        [Dict{String,Any}()], Dict{String,Any}()
    )
    
    lc = LightCurve(times, dt, counts, nothing, nothing, 
                    EventProperty{Float64}[], metadata, :poisson)
    flux_iterable = [counts]
    
    result_lightcurve = avg_pds_from_lightcurve(lc, silent=true)
    result_iterable = avg_pds_from_iterable(flux_iterable, dt, silent=true)
    
    @test result_lightcurve.power ≈ result_iterable.power
    @test result_lightcurve.freq ≈ result_iterable.freq
end

# Test LightCurve PDS function with different norms
let
    dt = 1.0
    N = 100
    times = collect(0:dt:(N-1)*dt)
    counts = rand(Poisson(100), N)
    
    metadata = LightCurveMetadata(
        "TEST", "INST", "TEST_OBJ", 58000.0, 
        (times[1], times[end]), dt,
        [Dict{String,Any}()], Dict{String,Any}()
    )
    
    lc = LightCurve(times, dt, counts, nothing, nothing, 
                    EventProperty{Float64}[], metadata, :poisson)
    flux_iterable = [counts]
    
    for norm in ["frac", "abs", "none", "leahy"]
        for use_common_mean in [true, false]
            result_lightcurve = avg_pds_from_lightcurve(lc, norm=norm, use_common_mean=use_common_mean, silent=true)
            result_iterable = avg_pds_from_iterable(flux_iterable, dt, norm=norm, use_common_mean=use_common_mean, silent=true)
            
            @test result_lightcurve.power ≈ result_iterable.power
            @test result_lightcurve.freq ≈ result_iterable.freq
        end
    end
end

# Test LightCurve cross spectrum function equivalence
let
    dt = 1.0
    N = 100
    times = collect(0:dt:(N-1)*dt)
    counts1 = rand(Poisson(100), N)
    counts2 = rand(Poisson(100), N)
    
    metadata = LightCurveMetadata(
        "TEST", "INST", "TEST_OBJ", 58000.0, 
        (times[1], times[end]), dt,
        [Dict{String,Any}()], Dict{String,Any}()
    )
    
    lc1 = LightCurve(times, dt, counts1, nothing, nothing, 
                     EventProperty{Float64}[], metadata, :poisson)
    lc2 = LightCurve(times, dt, counts2, nothing, nothing, 
                     EventProperty{Float64}[], metadata, :poisson)
    
    flux_iterable1 = [counts1]
    flux_iterable2 = [counts2]
    
    result_lightcurves = avg_cs_from_lightcurves(lc1, lc2, silent=true)
    result_iterables = avg_cs_from_iterables(flux_iterable1, flux_iterable2, dt, silent=true)
    
    @test result_lightcurves.power ≈ result_iterables.power
    @test result_lightcurves.freq ≈ result_iterables.freq
end

# Test LightCurve cross spectrum function with different parameters
let
    dt = 1.0
    N = 100
    times = collect(0:dt:(N-1)*dt)
    counts1 = rand(Poisson(100), N)
    counts2 = rand(Poisson(100), N)
    
    metadata = LightCurveMetadata(
        "TEST", "INST", "TEST_OBJ", 58000.0, 
        (times[1], times[end]), dt,
        [Dict{String,Any}()], Dict{String,Any}()
    )
    
    lc1 = LightCurve(times, dt, counts1, nothing, nothing, 
                     EventProperty{Float64}[], metadata, :poisson)
    lc2 = LightCurve(times, dt, counts2, nothing, nothing, 
                     EventProperty{Float64}[], metadata, :poisson)
    
    flux_iterable1 = [counts1]
    flux_iterable2 = [counts2]
    
    for norm in ["frac", "abs", "none", "leahy"]
        for use_common_mean in [true, false]
            for fullspec in [true, false]
                for return_auxil in [true, false]
                    result_lightcurves = avg_cs_from_lightcurves(
                        lc1, lc2, norm=norm, use_common_mean=use_common_mean,
                        fullspec=fullspec, return_auxil=return_auxil, silent=true
                    )
                    result_iterables = avg_cs_from_iterables(
                        flux_iterable1, flux_iterable2, dt,
                        norm=norm, use_common_mean=use_common_mean,
                        fullspec=fullspec, return_auxil=return_auxil, silent=true
                    )
                    
                    @test result_lightcurves.power ≈ result_iterables.power
                    @test result_lightcurves.freq ≈ result_iterables.freq
                    
                    if return_auxil
                        @test result_lightcurves.pds1 ≈ result_iterables.pds1
                        @test result_lightcurves.pds2 ≈ result_iterables.pds2
                    end
                end
            end
        end
    end
end

# Test EventList with no GTI (should create default GTI)
let
    data = setup_fourier_test_data()
    eventlist_no_gti = EventList(data.times, nothing, FITSMetadata(
        "[test]", 1, nothing, Dict{String,Vector}(), 
        Dict{String,Any}(), nothing, nothing
    ))
    
    expected_gti = Float64[minimum(data.times) maximum(data.times)]
    
    result_eventlist = avg_pds_from_eventlist(eventlist_no_gti, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    result_times = avg_pds_from_events(data.times, expected_gti, FOURIER_SEGMENT_SIZE, data.dt, silent=true)
    
    @test result_eventlist.power ≈ result_times.power
    @test result_eventlist.freq ≈ result_times.freq
end

# Test LightCurve with vector dt (should use first element)
let
    dt_vec = [1.0, 1.0, 1.0]
    N = 100
    times = collect(0:dt_vec[1]:(N-1)*dt_vec[1])
    counts = rand(Poisson(100), N)
    metadata = LightCurveMetadata(
        "TEST", "INST", "TEST_OBJ", 58000.0, 
        (times[1], times[end]), dt_vec[1],
        [Dict{String,Any}()], Dict{String,Any}()
    )
    lc_vec_dt = LightCurve(times, dt_vec, counts, nothing, nothing, 
                          EventProperty{Float64}[], metadata, :poisson)
    flux_iterable = [counts]
    result_lightcurve = avg_pds_from_lightcurve(lc_vec_dt, silent=true)
    result_iterable = avg_pds_from_iterable(flux_iterable, dt_vec[1], silent=true)
    @test result_lightcurve.power ≈ result_iterable.power
    @test result_lightcurve.freq ≈ result_iterable.freq
end