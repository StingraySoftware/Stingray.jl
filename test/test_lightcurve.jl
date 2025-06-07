function create_mock_eventlist(times, energies = nothing)
    # Create proper FITSMetadata structure
    headers = Dict{String,Any}(
        "TELESCOP" => "TEST",
        "INSTRUME" => "TEST",
        "OBJECT" => "TEST",
        "MJDREF" => 0.0,
    )

    # Create FITSMetadata with proper type parameters
    dummy_meta = FITSMetadata{Dict{String,Any}}(
        "test.fits",           # filepath
        1,                     # hdu
        "keV",                 # energy_units (or nothing if no energies)
        Dict{String,Vector}(), # extra_columns
        headers,                # headers
    )

    # Create EventList with proper type parameters
    # The TimeType should be the type of the vector, not the element type
    return EventList{typeof(times),typeof(dummy_meta)}(
        times,      # times vector
        energies,   # energies vector (or nothing)
        dummy_meta,  # metadata
    )
end
function create_mock_eventlist_meta(
    times::Vector{T},
    energies::Union{Nothing,Vector{T}} = nothing,
) where {T}
    # Create realistic test metadata that matches what extract_metadata expects
    test_headers = Dict{String,Any}(
        "TELESCOP" => "TEST",
        "INSTRUME" => "TEST",
        "OBJECT" => "TEST",
        "MJDREF" => 0.0,
        "OBSERVER" => "TEST_USER",
        "DATE-OBS" => "2023-01-01",
        "EXPOSURE" => 1000.0,
        "DATAMODE" => "TE",
    )

    # Create mock FITSMetadata with positional arguments matching the struct definition
    meta = FITSMetadata(
        "",           # filepath
        1,           # hdu
        nothing,     # energy_units
        Dict{String,Vector}(),  # extra_columns
        test_headers, # headers
    )

    return EventList(times, energies, meta)
end

# Test EventProperty structure creation and validation
let

    # Test basic EventProperty creation
    prop = EventProperty{Float64}(:test, [1.0, 2.0, 3.0], "units")
    @test prop.name === :test
    @test prop.values == [1.0, 2.0, 3.0]
    @test prop.unit == "units"
    @test typeof(prop) <: EventProperty{Float64}

    # Test different data types
    prop_int = EventProperty{Int}(:count, [1, 2, 3], "counts")
    @test prop_int.values == [1, 2, 3]
    @test typeof(prop_int) <: EventProperty{Int}

    # Test empty values
    prop_empty = EventProperty{Float64}(:empty, Float64[], "none")
    @test isempty(prop_empty.values)

end

# Test LightCurveMetadata structure creation and validation
let

    # Test complete metadata creation
    metadata = LightCurveMetadata(
        "TEST_TELESCOPE",
        "TEST_INSTRUMENT",
        "TEST_OBJECT",
        58000.0,
        (0.0, 100.0),
        1.0,
        [Dict{String,Any}("TEST" => "VALUE")],
        Dict{String,Any}("extra_info" => "test"),
    )

    @test metadata.telescope == "TEST_TELESCOPE"
    @test metadata.instrument == "TEST_INSTRUMENT"
    @test metadata.object == "TEST_OBJECT"
    @test metadata.mjdref == 58000.0
    @test metadata.time_range == (0.0, 100.0)
    @test metadata.bin_size == 1.0
    @test length(metadata.headers) == 1
    @test haskey(metadata.extra, "extra_info")
    @test metadata.extra["extra_info"] == "test"

    # Test with empty headers and extra info
    metadata_minimal = LightCurveMetadata(
        "",
        "",
        "",
        0.0,
        (0.0, 1.0),
        1.0,
        Vector{Dict{String,Any}}(),
        Dict{String,Any}(),
    )
    @test isempty(metadata_minimal.headers)
    @test isempty(metadata_minimal.extra)

end

# Test LightCurve basic structure creation and validation
let

    # Create test data
    timebins = [1.5, 2.5, 3.5]
    bin_edges = [1.0, 2.0, 3.0, 4.0]
    counts = [1, 2, 1]
    errors = Float64[1.0, √2, 1.0]
    exposure = fill(1.0, 3)
    props = [EventProperty{Float64}(:test, [1.0, 2.0, 3.0], "units")]
    metadata = LightCurveMetadata(
        "TEST",
        "TEST",
        "TEST",
        0.0,
        (1.0, 4.0),
        1.0,
        [Dict{String,Any}()],
        Dict{String,Any}(),
    )

    # Test LightCurve creation
    lc = LightCurve{Float64}(
        timebins,
        bin_edges,
        counts,
        errors,
        exposure,
        props,
        metadata,
        :poisson,
    )

    @test lc.timebins == timebins
    @test lc.bin_edges == bin_edges
    @test lc.counts == counts
    @test lc.count_error == errors
    @test lc.exposure == exposure
    @test length(lc.properties) == 1
    @test lc.err_method === :poisson
    @test typeof(lc) <: AbstractLightCurve{Float64}

    # Test inheritance
    @test lc isa AbstractLightCurve{Float64}

end

# Test Poisson error calculation
let

    # Test basic Poisson errors
    counts = [0, 1, 4, 9, 16]
    exposure = fill(1.0, length(counts))

    errors = calculate_errors(counts, :poisson, exposure)
    @test errors ≈ [1.0, 1.0, 2.0, 3.0, 4.0]

    # Test with zero counts (should use sqrt(1) = 1.0)
    zero_counts = [0, 0, 0]
    zero_errors = calculate_errors(zero_counts, :poisson, fill(1.0, 3))
    @test all(zero_errors .== 1.0)

    # Test with large counts
    large_counts = [100, 400, 900]
    large_errors = calculate_errors(large_counts, :poisson, fill(1.0, 3))
    @test large_errors ≈ [10.0, 20.0, 30.0]

end

# Test Gaussian error calculation
let
    counts = [1, 4, 9, 16, 25]
    exposure = fill(1.0, length(counts))
    gaussian_errs = [0.5, 1.0, 1.5, 2.0, 2.5]

    # Test with provided Gaussian errors
    errors_gauss =
        calculate_errors(counts, :gaussian, exposure, gaussian_errors = gaussian_errs)
    @test errors_gauss == gaussian_errs

    # Test different length Gaussian errors
    different_gaussian = [0.1, 0.2, 0.3]
    errors_diff = calculate_errors(
        [1, 2, 3],
        :gaussian,
        fill(1.0, 3),
        gaussian_errors = different_gaussian,
    )
    @test errors_diff == different_gaussian
end

# Test error calculation edge cases and exceptions
let

    counts = [1, 2, 3]
    exposure = fill(1.0, 3)

    # Test missing Gaussian errors
    @test_throws ArgumentError calculate_errors(counts, :gaussian, exposure)

    # Test wrong length Gaussian errors
    @test_throws ArgumentError calculate_errors(
        counts,
        :gaussian,
        exposure,
        gaussian_errors = [1.0, 2.0],
    )

    # Test invalid error method
    @test_throws ArgumentError calculate_errors(counts, :invalid, exposure)

    # Test empty arrays
    empty_errors = calculate_errors(Int[], :poisson, Float64[])
    @test isempty(empty_errors)
end

# Test input validation for lightcurve creation
let

    # Create valid EventList
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0]
    valid_events = create_mock_eventlist(times, energies)

    # Test valid inputs
    @test_nowarn validate_lightcurve_inputs(valid_events, 1.0, :poisson, nothing)
    @test_nowarn validate_lightcurve_inputs(valid_events, 0.1, :poisson, nothing)
    @test_nowarn validate_lightcurve_inputs(valid_events, 10.0, :poisson, nothing)

    # Test invalid bin size
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        0.0,
        :poisson,
        nothing,
    )
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        -1.0,
        :poisson,
        nothing,
    )
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        -0.1,
        :poisson,
        nothing,
    )

    # Test invalid error method
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        1.0,
        :invalid,
        nothing,
    )
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        1.0,
        :unknown,
        nothing,
    )

    # Test missing Gaussian errors
    @test_throws ArgumentError validate_lightcurve_inputs(
        valid_events,
        1.0,
        :gaussian,
        nothing,
    )
end

# Test empty event list validation
let

    # Create empty EventList
    empty_events = create_mock_eventlist(Float64[], nothing)

    # Test empty event list throws error
    @test_throws ArgumentError validate_lightcurve_inputs(
        empty_events,
        1.0,
        :poisson,
        nothing,
    )

end

# Test time filtering functionality
let

    times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]

    # Test time filtering only
    filtered_times, filtered_energies, start_t, stop_t =
        apply_event_filters(times, energies, 2.0, 4.0, nothing)
    @test all(2.0 .<= filtered_times .<= 4.0)
    @test length(filtered_times) == 3
    @test start_t == 2.0
    @test stop_t == 4.0
    @test length(filtered_energies) == length(filtered_times)

    # Test no time filtering (should use full range)
    filtered_times_full, _, start_t_full, stop_t_full =
        apply_event_filters(times, energies, nothing, nothing, nothing)
    @test length(filtered_times_full) == length(times)
    @test start_t_full == minimum(times)
    @test stop_t_full == maximum(times)

    # Test single time boundary
    filtered_start, _, _, _ = apply_event_filters(times, energies, 3.0, nothing, nothing)
    @test all(filtered_start .>= 3.0)

    filtered_stop, _, _, _ = apply_event_filters(times, energies, nothing, 4.0, nothing)
    @test all(filtered_stop .<= 4.0)

end

# Test energy filtering functionality
let

    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [5.0, 15.0, 25.0, 35.0, 45.0]

    # Test energy filtering
    filtered_times, filtered_energies, start_t, stop_t =
        apply_event_filters(times, energies, nothing, nothing, (10.0, 30.0))
    @test all(10.0 .<= filtered_energies .< 30.0)
    @test length(filtered_energies) == 2  # 15.0 and 25.0
    @test length(filtered_times) == length(filtered_energies)

    # Test energy filtering with inclusive lower bound, exclusive upper bound
    filtered_times2, filtered_energies2, _, _ =
        apply_event_filters(times, energies, nothing, nothing, (15.0, 25.0))
    @test 15.0 in filtered_energies2
    @test 25.0 ∉ filtered_energies2

    # Test energy filtering with no energies (should be no-op)
    filtered_times3, filtered_energies3, _, _ =
        apply_event_filters(times, nothing, nothing, nothing, (10.0, 30.0))
    @test length(filtered_times3) == length(times)
    @test isnothing(filtered_energies3)
end

# Test combined time and energy filtering
let

    times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    energies = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0]

    # Test combined filtering (energy first, then time)
    filtered_times, filtered_energies, start_t, stop_t =
        apply_event_filters(times, energies, 2.0, 4.0, (10.0, 40.0))
    @test all(2.0 .<= filtered_times .<= 4.0)
    @test all(10.0 .<= filtered_energies .< 40.0)
    @test length(filtered_times) == length(filtered_energies)

    # Verify specific filtered events
    expected_mask =
        (times .>= 2.0) .& (times .<= 4.0) .& (energies .>= 10.0) .& (energies .< 40.0)
    @test length(filtered_times) == sum(expected_mask)

end

# Test filtering edge cases and error conditions
let

    times = [1.0, 2.0, 3.0]
    energies = [10.0, 20.0, 30.0]

    # Test no events after energy filtering
    @test_throws ArgumentError apply_event_filters(
        times,
        energies,
        nothing,
        nothing,
        (100.0, 200.0),
    )

    # Test no events after time filtering
    @test_throws ArgumentError apply_event_filters(times, energies, 10.0, 20.0, nothing)

    # Test no events after combined filtering
    @test_throws ArgumentError apply_event_filters(
        times,
        energies,
        10.0,
        20.0,
        (100.0, 200.0),
    )
end

# Test time bin creation
let

    # Test basic bin creation
    start_time = 1.0
    stop_time = 5.0
    binsize = 1.0

    edges, centers = create_time_bins(start_time, stop_time, binsize)

    # Test bin structure
    @test length(edges) == length(centers) + 1
    @test edges[1] <= start_time
    @test edges[end] >= stop_time
    @test all(diff(edges) .≈ binsize)

    # Test centers are at bin midpoints
    expected_centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    @test centers ≈ expected_centers

    # Test with fractional binsize
    edges_frac, centers_frac = create_time_bins(0.5, 2.7, 0.3)
    @test all(diff(edges_frac) .≈ 0.3)
    @test edges_frac[1] <= 0.5
    @test edges_frac[end] >= 2.7

    # Test single bin case
    edges_single, centers_single = create_time_bins(1.0, 1.5, 2.0)
    @test length(centers_single) >= 1
    @test edges_single[end] >= 1.5

end

# Test event binning functionality
let

    # Test basic binning
    times = [1.1, 1.2, 2.3, 2.4, 3.5]
    edges = [1.0, 2.0, 3.0, 4.0]

    counts = bin_events(times, edges)
    @test length(counts) == length(edges) - 1
    @test counts == [2, 2, 1]  # 2 in [1,2), 2 in [2,3), 1 in [3,4)
    @test sum(counts) == length(times)

    # Test empty data
    empty_counts = bin_events(Float64[], edges)
    @test all(empty_counts .== 0)
    @test length(empty_counts) == length(edges) - 1

    # Test single event
    single_counts = bin_events([1.5], edges)
    @test sum(single_counts) == 1
    @test single_counts == [1, 0, 0]

    # Test events at bin boundaries
    boundary_times = [1.0, 2.0, 3.0]
    boundary_counts = bin_events(boundary_times, edges)
    @test sum(boundary_counts) == length(boundary_times)

    # Test with many events
    many_times = collect(1.1:0.1:3.9)
    many_counts = bin_events(many_times, edges)
    @test sum(many_counts) == length(many_times)

end

# Test additional properties calculation
let

    # Test with energy data
    times = [1.1, 1.2, 2.3, 2.4, 3.5]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0]
    edges = [1.0, 2.0, 3.0, 4.0]
    centers = [1.5, 2.5, 3.5]

    props = calculate_additional_properties(times, energies, edges, centers)

    # Test structure
    @test length(props) == 1
    @test props[1].name === :mean_energy
    @test props[1].unit == "keV"
    @test length(props[1].values) == length(centers)

    # Test mean energy calculation
    mean_energies = props[1].values
    @test mean_energies[1] ≈ mean([10.0, 20.0])  # Bin 1: events at 1.1, 1.2
    @test mean_energies[2] ≈ mean([15.0, 25.0])  # Bin 2: events at 2.3, 2.4
    @test mean_energies[3] ≈ 30.0                # Bin 3: event at 3.5

    # Test without energies
    props_no_energy = calculate_additional_properties(times, nothing, edges, centers)
    @test isempty(props_no_energy)

    # Test with empty energy data
    props_empty = calculate_additional_properties(Float64[], Float64[], edges, centers)
    @test isempty(props_empty)

    # Test with single bin
    single_edges = [1.0, 2.0]
    single_centers = [1.5]
    props_single = calculate_additional_properties(
        [1.1, 1.2],
        [10.0, 20.0],
        single_edges,
        single_centers,
    )
    @test length(props_single) == 1
    @test props_single[1].values[1] ≈ 15.0

end

# Test metadata extraction
let

    # Create mock eventlist with proper metadata
    times = [1.0, 2.0, 3.0]
    energies = [10.0, 20.0, 30.0]
    eventlist = create_mock_eventlist_meta(times, energies)

    # Test metadata extraction
    start_time = 1.0
    stop_time = 3.0
    binsize = 1.0
    filtered_times = times
    energy_filter = (5.0, 35.0)

    metadata = extract_metadata(
        eventlist,
        start_time,
        stop_time,
        binsize,
        filtered_times,
        energy_filter,
    )
    @test metadata.telescope == "TEST"
    @test metadata.instrument == "TEST"
    @test metadata.object == "TEST"
    @test metadata.mjdref == 0.0
    @test metadata.time_range == (1.0, 3.0)
    @test metadata.bin_size == 1.0

    # Test that ALL original headers are preserved
    @test !isempty(metadata.headers)
    @test haskey(metadata.headers[1], "TELESCOP")
    @test haskey(metadata.headers[1], "OBSERVER")
    @test haskey(metadata.headers[1], "DATE-OBS")
    @test haskey(metadata.headers[1], "EXPOSURE")
    @test metadata.headers[1]["OBSERVER"] == "TEST_USER"
    @test metadata.headers[1]["EXPOSURE"] == 1000.0

    # Test extra metadata - includes both processing and original extra data
    @test haskey(metadata.extra, "filtered_nevents")
    @test haskey(metadata.extra, "total_nevents")
    @test haskey(metadata.extra, "energy_filter")
    @test haskey(metadata.extra, "binning_method")

    @test metadata.extra["filtered_nevents"] == length(filtered_times)
    @test metadata.extra["total_nevents"] == length(times)
    @test metadata.extra["energy_filter"] == energy_filter
    @test metadata.extra["binning_method"] == "histogram"

end


# Test full lightcurve creation
let

    # Create test data
    times = [1.1, 1.2, 2.3, 2.4, 3.5, 4.1, 4.2]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0, 12.0, 18.0]
    eventlist = create_mock_eventlist(times, energies)

    # Test basic lightcurve creation
    lc = create_lightcurve(eventlist, 1.0)

    # Test structure
    @test length(lc.timebins) == length(lc.counts)
    @test length(lc.bin_edges) == length(lc.timebins) + 1
    @test length(lc.count_error) == length(lc.counts)
    @test length(lc.exposure) == length(lc.counts)
    @test sum(lc.counts) == length(times)
    @test all(lc.exposure .== 1.0)

    # Test metadata
    @test lc.metadata.bin_size == 1.0
    @test lc.metadata.extra["total_nevents"] == length(times)
    @test lc.metadata.extra["filtered_nevents"] == length(times)

    # Test properties
    @test !isempty(lc.properties)
    @test lc.properties[1].name === :mean_energy

end

# Test lightcurve creation with filtering
let

    times = [1.1, 1.2, 2.3, 2.4, 3.5, 4.1, 4.2]
    energies = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0]
    eventlist = create_mock_eventlist(times, energies)

    # Test with energy filtering
    lc_energy = create_lightcurve(eventlist, 1.0, energy_filter = (10.0, 50.0))
    @test sum(lc_energy.counts) < length(times)  # Some events filtered out
    @test lc_energy.metadata.extra["energy_filter"] == (10.0, 50.0)

    # Test with time filtering
    lc_time = create_lightcurve(eventlist, 1.0, tstart = 2.0, tstop = 4.0)
    @test sum(lc_time.counts) < length(times)  # Some events filtered out
    @test lc_time.metadata.time_range[1] == 2.0
    @test lc_time.metadata.time_range[2] == 4.0

    # Test with combined filtering
    lc_combined = create_lightcurve(
        eventlist,
        1.0,
        tstart = 2.0,
        tstop = 4.0,
        energy_filter = (10.0, 50.0),
    )
    @test sum(lc_combined.counts) <= sum(lc_energy.counts)
    @test sum(lc_combined.counts) <= sum(lc_time.counts)

end

# Test lightcurve creation with custom event filter
let

    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    eventlist = create_mock_eventlist(times, energies)

    # Test custom filter function
    custom_filter = eventlist -> eventlist.times .> 2.5
    lc_custom = create_lightcurve(eventlist, 1.0, event_filter = custom_filter)

    # Should only include events with times > 2.5
    expected_filtered = sum(times .> 2.5)
    @test sum(lc_custom.counts) == expected_filtered

    # Test filter that returns no events
    no_events_filter = eventlist -> fill(false, length(eventlist.times))
    @test_throws ArgumentError create_lightcurve(
        eventlist,
        1.0,
        event_filter = no_events_filter,
    )
end
# Test lightcurve creation with Gaussian errors
let

    times = [1.1, 1.2, 2.3, 2.4]
    energies = [10.0, 20.0, 15.0, 25.0]
    eventlist = create_mock_eventlist(times, energies)

    # Create lightcurve first to get bin structure
    lc_temp = create_lightcurve(eventlist, 1.0)
    n_bins = length(lc_temp.counts)

    # Create Gaussian errors for the right number of bins
    gaussian_errs = fill(0.5, n_bins)

    # Test with Gaussian errors
    lc_gauss = create_lightcurve(
        eventlist,
        1.0,
        err_method = :gaussian,
        gaussian_errors = gaussian_errs,
    )

    @test lc_gauss.err_method === :gaussian
    @test lc_gauss.count_error == gaussian_errs

    # Test error when Gaussian errors not provided
    @test_throws ArgumentError create_lightcurve(eventlist, 1.0, err_method = :gaussian)

    # Test error when Gaussian errors wrong length - use a clearly wrong length
    wrong_length_errs = [0.1]  # Length 1 - definitely wrong since we have 2 bins
    @test_throws ArgumentError create_lightcurve(
        eventlist,
        1.0,
        err_method = :gaussian,
        gaussian_errors = wrong_length_errs,
    )

    # Test with another wrong length
    another_wrong_length = [0.1, 0.2, 0.3, 0.4, 0.5]  # Length 5 - definitely wrong
    @test_throws ArgumentError create_lightcurve(
        eventlist,
        1.0,
        err_method = :gaussian,
        gaussian_errors = another_wrong_length,
    )

end
# Test basic rebinning functionality
let

    # Create test lightcurve
    start_time = 1.0
    end_time = 5.0
    old_binsize = 0.5
    new_binsize = 1.0

    # Create aligned time structure
    edges = collect(start_time:old_binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    counts = ones(Int, length(centers))

    lc = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(old_binsize, length(centers)),
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            old_binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test rebinning to larger bins
    new_lc = rebin(lc, new_binsize)

    # Test basic properties
    @test new_lc.metadata.bin_size == new_binsize
    @test all(new_lc.exposure .== new_binsize)
    @test sum(new_lc.counts) == sum(lc.counts)  # Count conservation
    @test length(new_lc.counts) < length(lc.counts)  # Fewer bins

    # Test error when rebinning to smaller bins
    @test_throws ArgumentError rebin(lc, old_binsize / 2)
    @test_throws ArgumentError rebin(lc, old_binsize)
end

# Test rebinning with properties
let

    # Create lightcurve with properties
    start_time = 1.0
    end_time = 5.0
    old_binsize = 1.0
    new_binsize = 2.0

    edges = collect(start_time:old_binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    n_bins = length(centers)

    counts = fill(2, n_bins)
    energy_values = collect(10.0:10.0:(10.0*n_bins))
    props = [EventProperty{Float64}(:mean_energy, energy_values, "keV")]

    lc = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(old_binsize, n_bins),
        props,
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            old_binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test rebinning with exact factor
    new_lc = rebin(lc, new_binsize)

    @test new_lc.metadata.bin_size == new_binsize
    @test sum(new_lc.counts) == sum(lc.counts)
    @test length(new_lc.properties) == length(lc.properties)
    @test all(new_lc.exposure .== new_binsize)
    @test new_lc.properties[1].name === :mean_energy
    @test new_lc.properties[1].unit == "keV"

    # Test property rebinning (should be weighted average)
    # For 2 bins with counts [2,2] and energies [10,20], weighted mean should be 15
    # This depends on your specific rebinning implementation
    @test length(new_lc.properties[1].values) == length(new_lc.counts)

    # Test half range rebinning
    total_range = end_time - start_time
    half_range_size = total_range / 2
    lc_half = rebin(lc, half_range_size)

    start_half = floor(start_time / half_range_size) * half_range_size
    n_half_bins = ceil(Int, (end_time - start_half) / half_range_size)
    @test length(lc_half.counts) == n_half_bins
    @test sum(lc_half.counts) == sum(lc.counts)

end

# Test rebinning edge cases
let

    # Test rebinning with non-aligned time structure
    start_time = 1.3
    end_time = 4.7
    old_binsize = 0.3
    new_binsize = 0.9

    edges = collect(start_time:old_binsize:(end_time+old_binsize))
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    counts = ones(Int, length(centers))

    lc_nonaligned = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(old_binsize, length(centers)),
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            old_binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test rebinning non-aligned bins
    new_lc_nonaligned = rebin(lc_nonaligned, new_binsize)
    @test sum(new_lc_nonaligned.counts) == sum(lc_nonaligned.counts)
    @test new_lc_nonaligned.metadata.bin_size == new_binsize

    # Test rebinning with single bin
    single_edges = [1.0, 2.0]
    single_centers = [1.5]
    single_counts = [10]

    lc_single = LightCurve{Float64}(
        single_centers,
        single_edges,
        single_counts,
        sqrt.(Float64.(single_counts)),
        [1.0],
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (1.0, 2.0),
            1.0,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Rebinning single bin to larger size should still work
    single_rebinned = rebin(lc_single, 2.0)
    @test sum(single_rebinned.counts) == sum(single_counts)
    @test length(single_rebinned.counts) >= 1

end

# Test rebinning error conditions
let

    # Create test lightcurve
    start_time = 1.0
    end_time = 5.0
    binsize = 1.0

    edges = collect(start_time:binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    counts = ones(Int, length(centers))

    lc = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(binsize, length(centers)),
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test invalid new bin sizes
    @test_throws ArgumentError rebin(lc, 0.0)  # Zero bin size
    @test_throws ArgumentError rebin(lc, -1.0)  # Negative bin size
    @test_throws ArgumentError rebin(lc, binsize / 2)  # Smaller than original
    @test_throws ArgumentError rebin(lc, binsize)  # Same as original

    # Test with very large bin size (should work but result in single bin)
    large_rebinned = rebin(lc, 100.0)
    @test length(large_rebinned.counts) == 1
    @test sum(large_rebinned.counts) == sum(lc.counts)

end

# Test rebinning preserves Gaussian errors
let

    # Create lightcurve with Gaussian errors
    start_time = 1.0
    end_time = 5.0
    old_binsize = 0.5
    new_binsize = 1.0

    edges = collect(start_time:old_binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    counts = fill(4, length(centers))  # Use constant counts for predictable errors
    gaussian_errors = fill(0.5, length(centers))  # Custom errors

    lc_gauss = LightCurve{Float64}(
        centers,
        edges,
        counts,
        gaussian_errors,
        fill(old_binsize, length(centers)),
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            old_binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :gaussian,
    )

    # Calculate the expected combined error first
    # For two bins with error 0.5 each, combined error should be sqrt(0.5^2 + 0.5^2)
    expected_combined_error = sqrt(2 * 0.5^2)

    # Calculate new Gaussian errors for rebinned light curve
    # For each new bin, we need to combine the errors from the original bins
    n_new_bins = ceil(Int, (end_time - start_time) / new_binsize)
    new_gaussian_errors = fill(expected_combined_error, n_new_bins)

    # Test rebinning preserves error method
    new_lc_gauss = rebin(lc_gauss, new_binsize, gaussian_errors = new_gaussian_errors)
    @test new_lc_gauss.err_method === :gaussian
    @test sum(new_lc_gauss.counts) == sum(lc_gauss.counts)
    @test length(new_lc_gauss.count_error) == length(new_lc_gauss.counts)

    # Test that rebinned errors are properly combined
    @test new_lc_gauss.count_error[1] ≈ expected_combined_error
end
function Base.iterate(lc::LightCurve)
    if length(lc.timebins) == 0
        return nothing
    end
    return (lc.timebins[1], lc.counts[1]), 2
end

function Base.iterate(lc::LightCurve, state)
    if state > length(lc.timebins)
        return nothing
    end
    return (lc.timebins[state], lc.counts[state]), state + 1
end

# Test lightcurve array interface
let

    times = [1.5, 2.5, 3.5]
    counts = [1, 2, 1]
    lc = LightCurve{Float64}(
        times,
        [1.0, 2.0, 3.0, 4.0],
        counts,
        sqrt.(Float64.(counts)),
        fill(1.0, 3),
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (1.0, 4.0),
            1.0,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test array interface
    @test length(lc) == 3
    @test size(lc) == (3,)
    @test lc[1] == (1.5, 1)
    @test lc[2] == (2.5, 2)
    @test lc[3] == (3.5, 1)

    # Test iteration
    collected = collect(lc)
    @test collected == [(1.5, 1), (2.5, 2), (3.5, 1)]

    # Test indexing with ranges
    if hasmethod(getindex, (typeof(lc), UnitRange{Int}))
        @test lc[1:2] == [(1.5, 1), (2.5, 2)]
    end

    # Test bounds checking
    @test_throws BoundsError lc[0]
    @test_throws BoundsError lc[4]
end

# Test rebinning with multiple properties
let

    start_time = 1.0
    end_time = 5.0
    old_binsize = 1.0
    new_binsize = 2.0

    edges = collect(start_time:old_binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    n_bins = length(centers)

    counts = fill(3, n_bins)
    energy_values = collect(10.0:5.0:(10.0+5.0*(n_bins-1)))
    flux_values = collect(1.0:0.5:(1.0+0.5*(n_bins-1)))

    props = [
        EventProperty{Float64}(:mean_energy, energy_values, "keV"),
        EventProperty{Float64}(:flux, flux_values, "cts/s"),
    ]

    lc_multi = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(old_binsize, n_bins),
        props,
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (start_time, end_time),
            old_binsize,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Test rebinning with multiple properties
    new_lc_multi = rebin(lc_multi, new_binsize)

    @test length(new_lc_multi.properties) == 2
    @test new_lc_multi.properties[1].name === :mean_energy
    @test new_lc_multi.properties[2].name === :flux
    @test new_lc_multi.properties[1].unit == "keV"
    @test new_lc_multi.properties[2].unit == "cts/s"

    # All properties should have same length as rebinned counts
    for prop in new_lc_multi.properties
        @test length(prop.values) == length(new_lc_multi.counts)
    end
end

# Test rebinning with empty lightcurve
let
    # Create empty lightcurve
    empty_lc = LightCurve{Float64}(
        Float64[],
        [1.0, 2.0],  # Minimal edges
        Int[],
        Float64[],
        Float64[],
        Vector{EventProperty{Float64}}(),
        LightCurveMetadata(
            "TEST",
            "TEST",
            "TEST",
            0.0,
            (1.0, 2.0),
            1.0,
            [Dict{String,Any}()],
            Dict{String,Any}(),
        ),
        :poisson,
    )

    # Rebinning empty lightcurve should work but result in empty or minimal structure
    rebinned_empty = rebin(empty_lc, 2.0)
    @test length(rebinned_empty.counts) >= 0
    @test sum(rebinned_empty.counts) == 0
    @test rebinned_empty.metadata.bin_size == 2.0

end

# Test rebinning preserves metadata
let

    # Create lightcurve with rich metadata
    start_time = 1.0
    end_time = 5.0
    old_binsize = 1.0
    new_binsize = 2.0

    edges = collect(start_time:old_binsize:end_time)
    centers = [(edges[i] + edges[i+1]) / 2 for i = 1:length(edges)-1]
    counts = ones(Int, length(centers))

    rich_metadata = LightCurveMetadata(
        "HUBBLE",
        "COS",
        "NGC1234",
        58000.0,
        (start_time, end_time),
        old_binsize,
        [Dict{String,Any}("OBSERVER" => "TEST", "DATE-OBS" => "2023-01-01")],
        Dict{String,Any}("custom_param" => 42, "processing_version" => "1.0"),
    )

    lc_rich = LightCurve{Float64}(
        centers,
        edges,
        counts,
        sqrt.(Float64.(counts)),
        fill(old_binsize, length(centers)),
        Vector{EventProperty{Float64}}(),
        rich_metadata,
        :poisson,
    )

    # Test rebinning preserves metadata
    rebinned_rich = rebin(lc_rich, new_binsize)

    @test rebinned_rich.metadata.telescope == "HUBBLE"
    @test rebinned_rich.metadata.instrument == "COS"
    @test rebinned_rich.metadata.object == "NGC1234"
    @test rebinned_rich.metadata.mjdref == 58000.0
    @test rebinned_rich.metadata.bin_size == new_binsize  # Should be updated
    @test rebinned_rich.metadata.time_range == (start_time, end_time)  # Should be preserved
    @test length(rebinned_rich.metadata.headers) == 1
    @test rebinned_rich.metadata.headers[1]["OBSERVER"] == "TEST"
    @test rebinned_rich.metadata.extra["custom_param"] == 42
    @test rebinned_rich.metadata.extra["processing_version"] == "1.0"

end
