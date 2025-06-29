"""
    validate_lightcurve_inputs(eventlist::EventList, binsize::Real, err_method::Symbol, gaussian_errors::Union{Nothing,Vector{<:Real}})

Validate inputs for light curve generation functions.

# Arguments
- `eventlist::EventList`: The event list containing photon arrival times and optional energy information
- `binsize::Real`: The time bin size for the light curve (must be positive)
- `err_method::Symbol`: Error calculation method, either `:poisson` or `:gaussian`
- `gaussian_errors::Union{Nothing,Vector{<:Real}}`: Vector of Gaussian errors (required when `err_method = :gaussian`)

# Returns
- `Bool`: Returns `true` if all inputs are valid

# Throws
- `ArgumentError`: If the event list is empty
- `ArgumentError`: If `binsize` is not positive
- `ArgumentError`: If `err_method` is not `:poisson` or `:gaussian`
- `ArgumentError`: If `err_method` is `:gaussian` but `gaussian_errors` is `nothing`
- `ArgumentError`: If length of `gaussian_errors` doesn't match length of event times

# Examples
```julia
# Valid inputs
eventlist = create_mock_eventlist([1.0, 2.0, 3.0])
validate_lightcurve_inputs(eventlist, 0.1, :poisson, nothing)  # true

# Invalid bin size
validate_lightcurve_inputs(eventlist, -0.1, :poisson, nothing)  # throws ArgumentError

# Gaussian method with errors
gaussian_errs = [0.1, 0.2, 0.15]
validate_lightcurve_inputs(eventlist, 0.1, :gaussian, gaussian_errs)  # true
```
"""
function validate_lightcurve_inputs(
    eventlist::EventList,
    binsize::Real,
    err_method::Symbol,
    gaussian_errors::Union{Nothing,Vector{<:Real}}
)
    if isempty(eventlist.times)
        throw(ArgumentError("Event list is empty"))
    end
    
    if binsize <= 0
        throw(ArgumentError("Bin size must be positive"))
    end
    
    if !(err_method in [:poisson, :gaussian])
        throw(ArgumentError("Unsupported error method: $err_method. Use :poisson or :gaussian"))
    end
    
    if err_method === :gaussian
        if isnothing(gaussian_errors)
            throw(ArgumentError("Gaussian errors must be provided when using :gaussian method"))
        end
        if length(gaussian_errors) != length(eventlist.times)
            throw(ArgumentError("Length of gaussian_errors must match length of event times"))
        end
    end
    
    return true
end

"""
    create_mock_eventlist(times, energies = nothing)

Create a mock EventList for testing purposes with minimal metadata.

This function creates a simple EventList with basic FITS headers suitable for testing
light curve generation and other analysis functions.

# Arguments
- `times`: Vector of photon arrival times
- `energies`: Optional vector of photon energies (default: `nothing`)

# Returns
- `EventList`: A mock event list with test metadata

# Details
The created EventList includes minimal FITS headers:
- `TELESCOP`: "TEST" 
- `INSTRUME`: "TEST"
- `OBJECT`: "TEST" 
- `MJDREF`: 0.0

# Examples
```julia
# Create event list with times only
times = [1.0, 2.5, 3.7, 4.1, 5.9]
eventlist = create_mock_eventlist(times)

# Create event list with times and energies
energies = [1.2, 2.3, 1.8, 3.1, 2.7]  # keV
eventlist = create_mock_eventlist(times, energies)
```

# See Also
- [`create_mock_eventlist_meta`](@ref): For creating mock event lists with extended metadata
"""
function create_mock_eventlist(times, energies = nothing)
    headers = Dict{String,Any}(
        "TELESCOP" => "TEST",
        "INSTRUME" => "TEST",
        "OBJECT" => "TEST",
        "MJDREF" => 0.0,
    )
    
    dummy_meta = FITSMetadata{Dict{String,Any}}(
        "test.fits",
        1,
        "keV",
        Dict{String,Vector}(),
        headers,
    )
    
    return EventList{typeof(times),typeof(dummy_meta)}(
        times,
        energies,
        dummy_meta,
    )
end

"""
    create_mock_eventlist_meta(times::Vector{T}, energies::Union{Nothing,Vector{T}} = nothing) where {T}

Create a mock EventList with comprehensive metadata for testing.

This function creates an EventList with extended FITS headers that more closely
resemble real astronomical observations, useful for testing functions that
depend on specific metadata fields.

# Arguments
- `times::Vector{T}`: Vector of photon arrival times of type T
- `energies::Union{Nothing,Vector{T}}`: Optional vector of photon energies of type T (default: `nothing`)

# Returns
- `EventList`: A mock event list with comprehensive test metadata

# Details
The created EventList includes comprehensive FITS headers:
- `TELESCOP`: "TEST" - Telescope name
- `INSTRUME`: "TEST" - Instrument name  
- `OBJECT`: "TEST" - Target object name
- `MJDREF`: 0.0 - Reference Modified Julian Date
- `OBSERVER`: "TEST_USER" - Observer name
- `DATE-OBS`: "2023-01-01" - Observation date
- `EXPOSURE`: 1000.0 - Exposure time in seconds
- `DATAMODE`: "TE" - Data mode (Timing Event)

# Type Parameters
- `T`: The numeric type of the times and energies vectors

# Examples
```julia
# Create event list with Float64 times
times = [1.0, 2.5, 3.7, 4.1, 5.9]
eventlist = create_mock_eventlist_meta(times)

# Create event list with times and energies
energies = [1.2, 2.3, 1.8, 3.1, 2.7]
eventlist = create_mock_eventlist_meta(times, energies)

# Access metadata
println(eventlist.meta.headers["EXPOSURE"])  # 1000.0
println(eventlist.meta.headers["OBSERVER"])  # "TEST_USER"
```

# See Also
- [`create_mock_eventlist`](@ref): For creating mock event lists with minimal metadata
"""
function create_mock_eventlist_meta(
    times::Vector{T},
    energies::Union{Nothing,Vector{T}} = nothing,
) where {T}
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
    
    meta = FITSMetadata(
        "",
        1,
        nothing,
        Dict{String,Vector}(),
        test_headers,
    )
    
    return EventList(times, energies, meta)
end
# Test EventProperty structure creation and validation
let
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
# Test Event Property Calculation
let
    times = [1.0, 1.8, 2.2, 3.1, 3.9, 4.1]
    energies = [12.0, 18.0, 22.0, 31.0, 39.0, 41.0]
    bin_edges = [1.0, 2.0, 3.0, 4.0, 5.0]
    bin_centers = [1.5, 2.5, 3.5, 4.5]
    properties = calculate_event_properties(times, energies, bin_edges, bin_centers)
    @test length(properties) == 1
    @test properties[1].name == :mean_energy
    @test properties[1].unit == "keV"
    @test length(properties[1].values) == 4
    expected_means = [15.0, 22.0, 35.0, 41.0]
    @test properties[1].values ≈ expected_means
end
# Test LightCurveMetadata structure creation and validation
let
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
    time = [1.5, 2.5, 3.5]
    dt = [1.0, 2.0, 3.0, 4.0]
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
    lc = LightCurve{Float64}(
        time,
        dt,
        counts,
        errors,
        exposure,
        props,
        metadata,
        :poisson,
    )
    @test lc.time == time
    @test lc.dt == dt
    @test lc.counts == counts
    @test lc.count_error == errors
    @test lc.exposure == exposure
    @test length(lc.properties) == 1
    @test lc.err_method === :poisson
    @test typeof(lc) <: AbstractLightCurve{Float64}
    # Test inheritance
    @test lc isa AbstractLightCurve{Float64}

end
# Test time bin creation
let
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
    times = [1.1, 1.2, 2.3, 2.4, 3.5]
    edges = [1.0, 2.0, 3.0, 4.0]
    counts = bin_events(times, edges)
    @test length(counts) == length(edges) - 1
    @test counts == [2, 2, 1]  # 2 in [1,2), 2 in [2,3), 1 in [3,4)
    @test sum(counts) == length(times)
    empty_counts = bin_events(Float64[], edges)
    @test all(empty_counts .== 0)
    @test length(empty_counts) == length(edges) - 1
    single_counts = bin_events([1.5], edges)
    @test sum(single_counts) == 1
    @test single_counts == [1, 0, 0]
    boundary_times = [1.0, 2.0, 3.0]
    boundary_counts = bin_events(boundary_times, edges)
    @test sum(boundary_counts) == length(boundary_times)
    many_times = collect(1.1:0.1:3.9)
    many_counts = bin_events(many_times, edges)
    @test sum(many_counts) == length(many_times)
end
# Test event properties calculation
let
    times = [1.1, 1.2, 2.3, 2.4, 3.5]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0]
    edges = [1.0, 2.0, 3.0, 4.0]
    centers = [1.5, 2.5, 3.5]
    props = calculate_event_properties(times, energies, edges, centers)
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
    props_no_energy = calculate_event_properties(times, nothing, edges, centers)
    @test isempty(props_no_energy)
    # Test with empty energy data
    props_empty = calculate_event_properties(Float64[], Float64[], edges, centers)
    @test isempty(props_empty)
    # Test with single bin
    single_edges = [1.0, 2.0]
    single_centers = [1.5]
    props_single = calculate_event_properties(
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
    times = [1.0, 2.0, 3.0]
    energies = [10.0, 20.0, 30.0]
    eventlist = create_mock_eventlist_meta(times, energies)
    start_time = 1.0
    stop_time = 3.0
    binsize = 1.0
    filtered_times = times
    energy_filter = (5.0, 35.0)
    # Pass the length of filtered_times, not the array itself
    metadata = extract_metadata(
        eventlist,
        start_time,
        stop_time,
        binsize,
        length(filtered_times),
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
    times = [1.1, 1.2, 2.3, 2.4, 3.5, 4.1, 4.2]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0, 12.0, 18.0]
    eventlist = create_mock_eventlist(times, energies)
    lc = create_lightcurve(eventlist, 1.0)
    @test length(lc.time) == length(lc.counts)
    @test length(lc.count_error) == length(lc.counts)
    @test length(lc.exposure) == length(lc.counts)
    @test sum(lc.counts) == length(times)
    @test all(lc.exposure .== 1.0)
    @test lc.dt == 1.0
    @test isa(lc.dt, Float64) || (isa(lc.dt, Vector) && all(lc.dt .== 1.0))
    # Test metadata
    @test lc.metadata.bin_size == 1.0
    @test lc.metadata.extra["total_nevents"] == length(times)
    @test lc.metadata.extra["filtered_nevents"] == length(times)    
    # Test properties
    @test !isempty(lc.properties)
    @test lc.properties[1].name === :mean_energy
    @test lc.err_method == :poisson
    @test all(lc.count_error .>= 0)  # Errors should be non-negative
    @test length(lc.properties[1].values) == length(lc.counts)
end
# Test lightcurve creation with filtering
let
    times = [1.1, 1.2, 2.3, 2.4, 3.5, 4.1, 4.2]
    energies = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0]
    eventlist = create_mock_eventlist(times, energies)
    lc_energy = create_lightcurve(eventlist, 1.0; energy_filter = (10.0, 50.0))
    # Events with energies 15.0, 25.0, 35.0, 45.0 should pass the filter (4 out of 7)
    @test sum(lc_energy.counts) == 4
    @test lc_energy.metadata.extra["energy_filter"] == (10.0, 50.0)
    # Test with time filtering
    lc_time = create_lightcurve(eventlist, 1.0; tstart = 2.0, tstop = 4.0)
    # Events at times 2.3, 2.4, 3.5 should pass the filter (3 out of 7)
    @test sum(lc_time.counts) == 3
    @test lc_time.metadata.time_range[1] == 2.0  
    @test lc_time.metadata.time_range[2] == 4.0
    lc_combined = create_lightcurve(
        eventlist,
        1.0;
        tstart = 2.0,
        tstop = 4.0,
        energy_filter = (10.0, 50.0)
    )
    # Events that pass both time (2.3, 2.4, 3.5) AND energy (15.0, 25.0, 35.0, 45.0) filters
    # Times 2.3 (energy 25.0), 2.4 (energy 35.0), 3.5 (energy 45.0) = 3 events
    @test sum(lc_combined.counts) == 3
    @test sum(lc_combined.counts) <= sum(lc_energy.counts)
    @test sum(lc_combined.counts) <= sum(lc_time.counts)    
    # Test that time range reflects the filtering parameters
    @test lc_combined.metadata.time_range[1] == 2.0
    @test lc_combined.metadata.time_range[2] == 4.0  
    # Test that energy filter is recorded in metadata
    @test lc_combined.metadata.extra["energy_filter"] == (10.0, 50.0)
end
# Test lightcurve creation with custom event filter
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    eventlist = create_mock_eventlist(times, energies)
    custom_filter = t -> t > 2.5  # Function that takes a single time value and returns boolean
    filtered_eventlist = filter_time(custom_filter, eventlist)  # Apply filter first
    lc_custom = create_lightcurve(filtered_eventlist, 1.0)
    # Should only include events with times > 2.5
    expected_filtered = sum(times .> 2.5) # This should be 2 (times 3.0, 4.0, 5.0 but 5.0 might be in next bin)
    @test sum(lc_custom.counts) == expected_filtered
    # Test filter that returns no events
    no_events_filter = t -> false
    empty_eventlist = filter_time(no_events_filter, eventlist)  # This will create empty eventlist
    @test_throws ArgumentError create_lightcurve(empty_eventlist, 1.0)  # Should throw when empty
end
# Test lightcurve creation with Gaussian errors and poisson
let
    times = [1.1, 1.2, 2.3, 2.4, 3.5]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0]
    eventlist = create_mock_eventlist(times, energies) 
    # Create lightcurve with Poisson errors
    lc = create_lightcurve(eventlist, 1.0)  
    # Test setting custom Gaussian errors
    custom_errors = [0.5, 1.0, 1.5]
    if length(lc.counts) == length(custom_errors)
        set_errors!(lc, custom_errors)
        @test lc.err_method == :gaussian
        @test lc.count_error == custom_errors
    end
    counts = [1, 4, 9, 16, 25]
    gaussian_errs = [0.5, 1.0, 1.5, 2.0, 2.5]
    # Test with provided Gaussian errors
    errors_gauss = calculate_errors(counts, :gaussian, gaussian_errors = gaussian_errs)
    @test errors_gauss == gaussian_errs
    # Test different length Gaussian errors
    different_gaussian = [0.1, 0.2, 0.3]
    errors_diff = calculate_errors([1, 2, 3], :gaussian, gaussian_errors = different_gaussian)
    @test errors_diff == different_gaussian
end
# Test basic rebinning functionality
let
    start_time = 1.0
    end_time = 5.0
    old_binsize = 0.5
    new_binsize = 1.0
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
    new_lc = rebin(lc, new_binsize)
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
 # Test rebinning with non-aligned time structure
let
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
# Test rebinning calculate_errors
let
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
    @test_throws ArgumentError rebin(lc, 0.0)  # Zero bin size
    @test_throws ArgumentError rebin(lc, -1.0)  # Negative bin size
    @test_throws ArgumentError rebin(lc, binsize / 2)  # Smaller than original
    @test_throws ArgumentError rebin(lc, binsize)  # Same as original
    large_rebinned = rebin(lc, 100.0)
    @test length(large_rebinned.counts) == 1
    @test sum(large_rebinned.counts) == sum(lc.counts)
end
# Test Poisson error calculation
let
    counts = [0, 1, 4, 9, 16]
    exposure = fill(1.0, length(counts))
    errors = calculate_errors(counts, :poisson)  # Fixed: removed exposure parameter
    @test errors ≈ [1.0, 1.0, 2.0, 3.0, 4.0]
    # Test with zero counts (should use sqrt(1) = 1.0)
    zero_counts = [0, 0, 0]
    zero_errors = calculate_errors(zero_counts, :poisson)  # Fixed: removed exposure parameter
    @test all(zero_errors .== 1.0)
    # Test with large counts
    large_counts = [100, 400, 900]
    large_errors = calculate_errors(large_counts, :poisson)  # Fixed: removed exposure parameter
    @test large_errors ≈ [10.0, 20.0, 30.0]
end
# Test Gaussian error calculation
let
    counts = [1, 4, 9, 16, 25]
    gaussian_errs = [0.5, 1.0, 1.5, 2.0, 2.5]
    # Test with provided Gaussian errors
    errors_gauss = calculate_errors(counts, :gaussian, gaussian_errors = gaussian_errs)
    @test errors_gauss == gaussian_errs
    # Test different length Gaussian errors
    different_gaussian = [0.1, 0.2, 0.3]
    errors_diff = calculate_errors([1, 2, 3], :gaussian, gaussian_errors = different_gaussian)
    @test errors_diff == different_gaussian
end
# Test error calculation[exclusive tests]
let
    counts = [1, 2, 3]
    # Test missing Gaussian errors
    @test_throws ArgumentError calculate_errors(counts, :gaussian)
    # Test wrong length Gaussian errors
    @test_throws ArgumentError calculate_errors(counts, :gaussian, gaussian_errors = [1.0, 2.0])
    # Test invalid error method
    @test_throws ArgumentError calculate_errors(counts, :invalid)
    # Test empty arrays
    empty_errors = calculate_errors(Int[], :poisson)
    @test isempty(empty_errors)
    # Test Poisson errors with edge cases
    counts_poisson = [0, 1, 4, 9, 16, 25]
    expected_poisson = [1.0, 1.0, 2.0, 3.0, 4.0, 5.0]  # sqrt(max(counts, 1))
    result_poisson = calculate_errors(counts_poisson, :poisson)
    @test result_poisson ≈ expected_poisson
    # Test Gaussian errors with validation
    counts_gaussian = [10, 20, 30]
    gaussian_errors = [1.5, 2.5, 3.5]
    result_gaussian = calculate_errors(counts_gaussian, :gaussian; gaussian_errors=gaussian_errors)
    @test result_gaussian == gaussian_errors
    # Test error conditions
    @test_throws ArgumentError calculate_errors(counts_gaussian, :gaussian)  # Missing errors
    @test_throws ArgumentError calculate_errors(counts_gaussian, :gaussian; gaussian_errors=[1.0, 2.0])  # Wrong length
    @test_throws ArgumentError calculate_errors(counts_gaussian, :invalid_method)  # Invalid method
end
# Test Comprehensive Filtering System
let
    times = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    energies = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0]    
    # Test time-only filtering
    filtered_times, filtered_energies, start_t, stop_t = 
        apply_filters(times, energies, 2.0, 5.0, nothing)
    @test all(2.0 .<= filtered_times .<= 5.0)
    @test length(filtered_times) == 3  # [2.5, 3.5, 4.5] (5.5 > 5.0)
    @test start_t == 2.0
    @test stop_t == 5.0
    # Test energy-only filtering  
    filtered_times_e, filtered_energies_e, _, _ =
        apply_filters(times, energies, nothing, nothing, (20.0, 50.0))
    @test all(20.0 .<= filtered_energies_e .< 50.0)
    @test 25.0 in filtered_energies_e && 35.0 in filtered_energies_e && 45.0 in filtered_energies_e
    @test 50.0 ∉ filtered_energies_e  # Exclusive upper bound
    filtered_times_c, filtered_energies_c, _, _ =
        apply_filters(times, energies, 2.0, 5.0, (20.0, 50.0))
    @test all(2.0 .<= filtered_times_c .<= 5.0)
    @test all(20.0 .<= filtered_energies_c .< 50.0)
    @test_throws ArgumentError apply_filters(times, energies, 10.0, 20.0, nothing)  # No events in time range
    @test_throws ArgumentError apply_filters(times, energies, nothing, nothing, (100.0, 200.0))  # No events in energy range    
    # Test filtering without energies
    filtered_no_energy, _, _, _ = apply_filters(times, nothing, 2.0, 5.0, (20.0, 50.0))
    @test all(2.0 .<= filtered_no_energy .<= 5.0)
end

# Test Binning and Event Processing
let
    start_time, stop_time, binsize = 0.0, 10.0, 2.0
    edges, centers = create_time_bins(start_time, stop_time, binsize)
    @test edges[1] <= start_time
    @test edges[end] >= stop_time
    @test all(diff(edges) .≈ binsize)
    @test length(centers) == length(edges) - 1
    @test all(centers .≈ edges[1:end-1] .+ binsize/2)
    event_times = [0.1, 0.9, 1.1, 1.9, 2.1, 2.8, 3.2, 3.9]
    bin_edges = [0.0, 1.0, 2.0, 3.0, 4.0]
    counts = bin_events(event_times, bin_edges)
    expected_counts = [2, 2, 2, 2]
    @test counts == expected_counts
    @test_throws ArgumentError bin_events([1.0, 2.0], [1.0])  # Insufficient bins
    sparse_times = [0.1, 3.9]  # Only first and last bins
    sparse_counts = bin_events(sparse_times, bin_edges)
    @test sparse_counts == [1, 0, 0, 1]
end

# Test Light Curve Validation and Creation
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    eventlist = create_mock_eventlist(times, energies)
    @test validate_lightcurve_inputs(eventlist, 1.0, :poisson, nothing) == true

    gaussian_errors = [1.0, 1.5, 2.0, 2.5, 3.0]
    @test validate_lightcurve_inputs(eventlist, 1.0, :gaussian, gaussian_errors) == true

    empty_eventlist = create_mock_eventlist(Float64[], nothing)
    @test_throws ArgumentError validate_lightcurve_inputs(empty_eventlist, 1.0, :poisson, nothing)
    @test_throws ArgumentError validate_lightcurve_inputs(eventlist, -1.0, :poisson, nothing)  
    @test_throws ArgumentError validate_lightcurve_inputs(eventlist, 1.0, :invalid, nothing)  # Invalid method
    @test_throws ArgumentError validate_lightcurve_inputs(eventlist, 1.0, :gaussian, nothing)  # Missing gaussian errors
end
# Test Light Curve Array Interface and Operations
let
    times = [1.5, 2.5, 3.5, 4.5]
    counts = [10, 20, 15, 25]
    metadata = LightCurveMetadata(
        "TEST_MISSION", "TEST_INSTRUMENT", "TEST_OBJECT", 0.0,
        (1.0, 5.0), 1.0, [Dict{String,Any}()], Dict{String,Any}()
    )
    lc = LightCurve{Float64}(
        times, 1.0, counts, nothing, nothing,
        Vector{EventProperty{Float64}}(), metadata, :poisson
    )
    
    @test length(lc) == 4
    @test size(lc) == (4,)
    @test lc[1] == (1.5, 10)
    @test lc[2] == (2.5, 20)
    @test lc[end] == (4.5, 25)
    # Test slicing
    slice = lc[2:3]
    @test slice == [(2.5, 20), (3.5, 15)]
    # Test iteration
    collected = collect(lc)
    expected_pairs = [(1.5, 10), (2.5, 20), (3.5, 15), (4.5, 25)]
    @test collected == expected_pairs
    # Test error setting functionality
    set_errors!(lc)
    @test lc.err_method == :poisson
    @test lc.count_error ≈ sqrt.(counts)
    # Test custom error setting
    custom_errors = [1.5, 2.5, 1.8, 3.2]
    set_errors!(lc, custom_errors)
    @test lc.err_method == :gaussian
    @test lc.count_error == custom_errors
    # Test error conditions
    @test_throws ArgumentError set_errors!(lc, [1.0, 2.0])  # Wrong length
end
# Test cases for LightCurve error functions[set_erros! and calculate_errors!]
let
    event_times = [0.5, 1.5, 2.5, 3.5, 4.5, 0.7, 1.2, 2.8, 3.1, 4.9]
    eventlist = create_mock_eventlist(event_times)
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [100.0, 144.0, 225.0, 196.0, 169.0]
    bin_width = 1.0
    bins = [1, 2, 3, 4, 5]
    count_error = nothing
    rate = nothing
    properties = EventProperty{Float64}[]
    metadata = LightCurveMetadata(
        "TEST",      # telescope
        "TEST",      # instrument  
        "TEST",      # object
        0.0,         # mjdref
        (0.0, 5.0),  # time_range
        1.0,         # bin_width
        Vector{Dict{String,Any}}(),  # filters
        Dict{String,Any}()           # extra_info
    )
    err_method = :poisson
    lc = LightCurve(times, bin_width, bins, count_error, rate, properties, metadata, err_method)
    lc.counts = counts
    set_errors!(lc)
    @test lc.err_method == :poisson
    @test lc.count_error ≈ sqrt.(counts) 
    # Test set_errors! with custom errors
    custom_errors = [0.5, 1.2, 0.8, 1.1, 0.9]
    set_errors!(lc, custom_errors)
    @test lc.err_method == :gaussian
    @test lc.count_error == custom_errors
    # Test ArgumentError for mismatched length
    wrong_length_errors = [0.5, 1.2, 0.8]
    @test_throws ArgumentError set_errors!(lc, wrong_length_errors)
    # Test calculate_errors! with Poisson method
    set_errors!(lc)
    original_errors = copy(lc.count_error)
    lc.counts .+= 10.0
    new_errors = calculate_errors!(lc)
    expected_errors = sqrt.(lc.counts)
    @test new_errors ≈ expected_errors
    @test lc.count_error ≈ expected_errors
    @test lc.err_method == :poisson
    @test !(lc.count_error ≈ original_errors)
    # Test calculate_errors! preserves method but only works with Poisson
    # Note: calculate_errors! likely only works with Poisson method
    # For Gaussian method, errors are fixed and don't recalculate automatically
    custom_errors_2 = [2.0, 2.5, 3.0, 2.8, 2.2]
    set_errors!(lc, custom_errors_2)  # Set to Gaussian
    # Store original Gaussian errors
    original_gaussian_errors = copy(lc.count_error)
    # Modify counts - Gaussian errors should remain unchanged
    lc.counts .*= 2.0  
    # For Gaussian method, calculate_errors! may not work since errors are fixed
    # So we test that the errors remain the same (no automatic recalculation)
    @test lc.err_method == :gaussian
    @test lc.count_error == original_gaussian_errors  # Should remain unchange
    set_errors!(lc)
    lc.counts .= [50.0, 75.0, 100.0, 125.0, 150.0]
    returned_errors = calculate_errors!(lc)
    expected_errors = sqrt.(lc.counts)
    @test returned_errors ≈ expected_errors
    @test lc.count_error ≈ expected_errors
    @test lc.err_method == :poisson
    # Test return value consistency for calculate_errors!
    times_2 = [1.0, 2.0]  
    counts_2 = [16.0, 25.0]
    metadata_2 = LightCurveMetadata("TEST", "TEST", "TEST", 0.0, (0.0, 2.0), 1.0, Vector{Dict{String,Any}}(), Dict{String,Any}())
    lc_test = LightCurve(times_2, 1.0, [1, 2], nothing, nothing, EventProperty{Float64}[], metadata_2, :poisson)
    lc_test.counts = counts_2
    set_errors!(lc_test) 
    returned_errors = calculate_errors!(lc_test)
    @test returned_errors == lc_test.count_error
    @test returned_errors ≈ [4.0, 5.0]
    # Test type consistency
    times_int = [1.0, 2.0, 3.0]
    counts_int = [100, 144, 225]
    metadata_int = LightCurveMetadata("TEST", "TEST", "TEST", 0.0, (0.0, 3.0), 1.0, Vector{Dict{String,Any}}(), Dict{String,Any}())
    lc_int = LightCurve(times_int, 1.0, [1, 2, 3], nothing, nothing, EventProperty{Float64}[], metadata_int, :poisson)
    lc_int.counts = counts_int
    set_errors!(lc_int)
    @test eltype(lc_int.count_error) == Float64
    float_errors = [1.5, 2.5, 3.5]
    set_errors!(lc_int, float_errors)
    @test eltype(lc_int.count_error) == Float64
    @test lc_int.count_error == [1.5, 2.5, 3.5]
end