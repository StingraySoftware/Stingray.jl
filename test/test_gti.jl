# Helper function to create mock EventList for testing
function create_test_eventlist(times::Vector{Float64}, energies::Union{Vector{Float64}, Nothing}=nothing)
    mock_headers = Dict{String,Any}()
    mock_metadata = FITSMetadata(
        "",  # filepath
        1,   # hdu
        "keV",  # energy_units
        Dict{String,Vector}(),  # extra_columns
        mock_headers,  # headers
        nothing,       # gti
        nothing        # gti_source
    )
    
    return EventList(times, energies, mock_metadata)
end

# Helper function to create mock LightCurve for testing
function create_test_lightcurve(times::Vector{Float64}, counts::Vector{Int}, dt::Float64=1.0)
    metadata = LightCurveMetadata(
        "", "", "", 0.0, 
        (minimum(times)-dt/2, maximum(times)+dt/2), 
        dt, 
        Vector{Dict{String,Any}}(),
        Dict{String,Any}()
    )
    
    return LightCurve(
        times, dt, counts, nothing, nothing, EventProperty{Float64}[], 
        metadata, :poisson
    )
end

# Helper function to create EventProperty
function create_event_property(name::String, values::Vector{Float64}, unit::String="")
    return EventProperty{Float64}(Symbol(name), values, unit)
end

# Helper functions for common test data
function get_basic_times()
    return [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
end

function get_basic_energies()
    return [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
end

function get_basic_counts()
    return [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
end

function get_simple_times()
    return [1.0, 2.0, 3.0, 4.0, 5.0]
end

function get_simple_energies()
    return [1.0, 2.0, 3.0, 4.0, 5.0]
end

function get_simple_counts()
    return [10, 20, 30, 40, 50]
end

function get_sparse_times()
    return [0.5, 1.5, 8.5, 9.5]
end

function get_sparse_energies()
    return [1.0, 2.0, 3.0, 4.0]
end

function get_bti_test_times()
    return [1.0, 2.0, 6.0, 7.0]
end

function get_bti_test_energies()
    return [1.0, 2.0, 3.0, 4.0]
end
# test_load_gtis
let
    fname = joinpath(@__DIR__ ,"data","monol_testA.evt")
    @test load_gtis(fname) == [8.0e7 8.0001025e7]
end

# get_gti_length
let
    @test get_total_gti_length([[0 5]; [6 7];]) == 6
end

# test_check_gtis_shape
let
    @test_throws ArgumentError check_gtis([[0 1 4]; [0 3 4];]) 
end

# test_check_gtis_values
let
    @test_throws ArgumentError check_gtis([[0 2]; [1 3];])

    @test_throws ArgumentError check_gtis([[1 0];])
end

# test_gti_mask1
let
    arr = [0, 1, 2, 3, 4, 5, 6]
    gti = [[0 2.1]; [3.9 5];]
    mask, new_gtis = create_gti_mask(arr, gti)
    # NOTE: the time bin has to be fully inside the GTI. That is why the
    # bin at times 0, 2, 4 and 5 are not in.
    @test mask == [0, 1, 0, 0, 0, 0, 0]
end

# test_gti_mask_minlen
let
    arr = [0, 1, 2, 3, 4, 5, 6]
    gti = [[0 2.1]; [3.9 5];]
    mask, new_gtis = create_gti_mask(arr, gti; min_length=2)
    # NOTE: the time bin has to be fully inside the GTI. That is why the
    # bin at times 0, 2, 4 and 5 are not in.
    @test mask == [0, 1, 0, 0, 0, 0, 0]
    @test new_gtis == [0 2.1]
end

# test_gti_mask_none_longer_than_minlen
let
    arr = [0, 1, 2, 3, 4, 5, 6]
    gti = [[0 2.1]; [3.9 5];]
    mask = Bool[]
    @test_logs (:warn,r"No GTIs longer than"
        ) mask, _ = create_gti_mask(arr, gti; min_length=10)
    @test all(iszero, mask)
end

# test_gti_mask_fails_empty_time
let
    arr = Float64[]
    gti = [[0 2.1]; [3.9 5];]
    @test_throws ArgumentError create_gti_mask(arr, gti)
end

# test_gti_from_condition1
let
    t = [0, 1, 2, 3, 4, 5, 6]
    condition = [true, true, false, false, true, false, false]
    gti = create_gti_from_condition(t, condition)
    @test gti == [[-0.5 1.5]; [3.5 4.5];]
end

# test_gti_from_condition2
let
    t = [0, 1, 2, 3, 4, 5, 6]
    condition = [true, true, true, true, false, true, false]
    gti = create_gti_from_condition(t, condition, safe_interval=[1, 1])
    @test gti == [0.5 2.5]
end

# test_gti_from_condition_fail
let
    t = [0, 1, 2, 3]
    condition = [true, true, true]
    @test_throws ArgumentError create_gti_from_condition(t, condition, safe_interval=[1, 1])
end

# test_intersectgti1
let
    gti1 = [1 4]
    gti2 = [2 5]
    newgti = operations_on_gtis([gti1, gti2], intersect)
    @test newgti == [2 4]
end

# test_intersectgti2
let
    gti1 = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
    gti2 = [[2 5]; [6 9]; [11.4 14];]
    newgti = operations_on_gtis([gti1, gti2], intersect)
    @test newgti == [[4.0 5.0]; [7.0 9.0]; [12.2 13.2];]
end

# test_intersectgti3
let
    gti1 = [[1 2]; [4 5]; [7 10];]
    newgti = operations_on_gtis([gti1], intersect)
    @test newgti == gti1
end

# test_union_gtis_nonoverlapping
let
    gti0 = [[0 1]; [2 3];]
    gti1 = [[10 11]; [12 13];]
    @test operations_on_gtis([gti0, gti1], union) == [[0 1]; [2 3]; [10 11]; [12 13];]
end

# test_union_gtis_overlapping
let
    gti0 = [[0 1]; [2 3]; [4 8];]
    gti1 = [[7 8]; [10 11]; [12 13];]
    @test operations_on_gtis([gti0, gti1], union) == [[0 1]; [2 3]; [4 8]; [10 11]; [12 13];]
end

# test_bti
let
    gti = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
    bti = get_btis(gti)
    @test bti == [[2 4]; [5 7]; [10 11]; [11.2 12.2];]
end

# test_bti_start_and_stop
let
    gti = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
    bti = get_btis(gti, 0, 14)
    @test bti == [[0 1]; [2 4]; [5 7]; [10 11]; [11.2 12.2]; [13.2 14];]
end

# test_bti_empty_valid
let
    gti = reshape(Float64[],0,2)
    bti = get_btis(gti, 0, 1)
    @test bti == [0 1]
end

# test_bti_fail
let
    gti = reshape(Float64[],0,2)
    @test_throws ArgumentError get_btis(gti)
end

# test_time_intervals_from_gtis1
let
    start_times, stop_times = time_intervals_from_gtis([[0 400]; [1022 1200];
                              [1210 1220];], 128)
    @test start_times == [0, 128, 256, 1022]
    @test stop_times == start_times .+ 128
end

# test_time_intervals_from_gtis_frac
let
    start_times, stop_times = time_intervals_from_gtis([[0 400]; [1022 1200];
                              [1210 1220];], 128, fraction_step=0.5)
    @test start_times == [0, 64, 128, 192, 256, 1022]
    @test stop_times == start_times .+ 128
end

# test_bin_intervals_from_gtis1
let
    times = range(0.5, 12.5)
    start_bins, stop_bins = bin_intervals_from_gtis([[0 5]; [6 8];], 2, times)

    @test start_bins == [0, 2, 6]
    @test stop_bins == [2, 4, 8]
end

# test_bin_intervals_from_gtis_2
let
    dt = 0.1
    tstart = 0
    tstop = 100
    times = range(tstart, tstop, step=dt)
    gti = [[tstart - dt/2 tstop - dt/2];]
    start_bins, stop_bins = bin_intervals_from_gtis(gti, 20, times)
    @test start_bins == [0, 200, 400, 600, 800]
end

# test_bin_intervals_from_gtis_frac
let
    times = range(0.5, 12.5)
    start_bins, stop_bins = bin_intervals_from_gtis([[0 5]; [6 8];], 2, times, fraction_step=0.5)

    @test start_bins == [0, 1, 2, 3, 6]
    @test stop_bins == [2, 3, 4, 5, 8]
end
# EventList GTI Tests - Basic functionality
let
    times = get_basic_times()
    energies = get_basic_energies()
    
    el = create_test_eventlist(times, energies)
    
    gtis = [2.0 6.0; 7.0 9.0]
    result = split_by_gtis(el, gtis)
    
    @test length(result) == 2
    @test result[1].times ≈ [2.5, 3.5, 4.5, 5.5]
    @test result[1].energies ≈ [3.0, 4.0, 5.0, 6.0]
    @test result[2].times ≈ [7.5, 8.5]
    @test result[2].energies ≈ [8.0, 9.0]
end

# EventList GTI Tests - No events in GTI
let
    times = get_sparse_times()
    energies = get_sparse_energies()
    
    el = create_test_eventlist(times, energies)
    gtis = [2.0 3.0; 4.0 5.0]
    
    result = split_by_gtis(el, gtis)
    @test length(result) == 0
end

# EventList GTI Tests - All events within single GTI
let
    times = get_simple_times()
    energies = get_simple_energies()
    
    el = create_test_eventlist(times, energies)
    gtis = [0.5 5.5]
    
    result = split_by_gtis(el, gtis)
    @test length(result) == 1
    @test result[1].times == times
    @test result[1].energies == energies
end

# LightCurve GTI Tests - Basic functionality
let
    time_bins = get_basic_times()
    counts = get_basic_counts()
    
    lc = create_test_lightcurve(time_bins, counts)
    
    gtis = [2.0 6.0; 7.0 9.0]
    result = apply_gtis(lc, gtis)
    
    @test length(result) == 2
    @test result[1].time ≈ [2.5, 3.5, 4.5, 5.5]
    @test result[1].counts ≈ [20, 25, 30, 35]
    @test result[2].time ≈ [7.5, 8.5]
    @test result[2].counts ≈ [45, 50]
end

# LightCurve GTI Tests - With exposure
let
    time_bins = get_simple_times()
    counts = get_simple_counts()
    exposure = [0.9, 1.0, 1.0, 0.8, 1.0]
    
    metadata = LightCurveMetadata(
        "", "", "", 0.0, (0.0, 6.0), 1.0, 
        Vector{Dict{String,Any}}(), Dict{String,Any}()
    )
    
    lc = LightCurve(
        time_bins, 1.0, counts, nothing, exposure, EventProperty{Float64}[], 
        metadata, :poisson
    )
    
    gtis = [1.5 3.5]
    result = apply_gtis(lc, gtis)
    
    @test length(result) == 1
    @test result[1].time ≈ [2.0, 3.0]
    @test result[1].counts ≈ [20, 30]
    @test result[1].exposure ≈ [1.0, 1.0]
end

# LightCurve GTI Tests - With properties
let
    time_bins = get_simple_times()
    counts = get_simple_counts()
    
    energy_prop = create_event_property("energy", [1.5, 2.5, 3.5, 4.5, 5.5], "keV")
    properties = [energy_prop]
    
    metadata = LightCurveMetadata(
        "", "", "", 0.0, (0.0, 6.0), 1.0, 
        Vector{Dict{String,Any}}(), Dict{String,Any}()
    )
    
    lc = LightCurve(
        time_bins, 1.0, counts, nothing, nothing, properties, 
        metadata, :poisson
    )
    
    gtis = [1.5 3.5]
    result = apply_gtis(lc, gtis)
    
    @test length(result) == 1
    @test result[1].properties[1].values ≈ [2.5, 3.5]
end

# Filtered LightCurve Tests - Basic filtering
let
    time_bins = get_simple_times()
    counts = get_simple_counts()
    
    lc = create_test_lightcurve(time_bins, counts)
    lc.metadata = LightCurveMetadata(
        "TEST", "INSTR", "OBJ", 58000.0, (0.0, 6.0), 1.0, 
        Vector{Dict{String,Any}}([Dict{String,Any}("OBSERVER" => "Test")]), 
        Dict{String,Any}("original_key" => "value")
    )
    
    mask = [false, true, true, false, true]
    filtered_lc = create_filtered_lightcurve(lc, mask, 1.5, 3.5, 1)
    
    @test filtered_lc.time ≈ [2.0, 3.0, 5.0]
    @test filtered_lc.counts ≈ [20, 30, 50]
    @test filtered_lc.dt == 1.0
    @test filtered_lc.metadata.extra["gti_applied"] == true
    @test filtered_lc.metadata.extra["gti_index"] == 1
    @test filtered_lc.metadata.extra["filtered_nbins"] == 3
    @test filtered_lc.metadata.extra["original_nbins"] == 5
end

# BTI Filling Tests - Basic BTI filling
let
    times = get_simple_times()
    energies = get_simple_energies()
    
    el = create_test_eventlist(times, energies)
    gtis = [1.0 2.5; 5.5 8.5]
    
    fill_bad_time_intervals!(el, gtis, dt=0.5, random_fill_threshold=5.0, 
                           rng=Random.MersenneTwister(123))
    
    @test length(el.times) > 5
    @test el.meta.headers["BTI_FILLED"] == true
    @test el.meta.headers["N_SYNTH_EVENTS"] > 0
    @test el.meta.headers["RAND_FILL_THRESH"] == 5.0
    @test el.meta.headers["BTI_FILL_DT"] == 0.5
    
    @test haskey(el.meta.extra_columns, "filled_bti_durations")
    @test length(el.meta.extra_columns["filled_bti_durations"]) > 0
    
    @test issorted(el.times)
end

# BTI Filling Tests - No BTI to fill
let
    times = [1.0, 2.0, 15.0, 16.0]
    energies = [1.0, 2.0, 3.0, 4.0]
    
    el = create_test_eventlist(times, energies)
    gtis = [1.0 2.5; 14.5 16.5]
    
    original_length = length(el.times)
    fill_bad_time_intervals!(el, gtis, dt=1.0, random_fill_threshold=10.0,
                           rng=Random.MersenneTwister(123))
    
    @test length(el.times) == original_length
    @test !haskey(el.meta.headers, "BTI_FILLED")
end

# BTI Filling Tests - Without energies
let
    times = get_bti_test_times()
    
    el = create_test_eventlist(times, nothing)
    gtis = [1.0 2.5; 5.5 7.5]
    
    fill_bad_time_intervals!(el, gtis, dt=1.0, random_fill_threshold=5.0,
                           rng=Random.MersenneTwister(123))
    
    @test el.energies === nothing
    @test length(el.times) > 4
end

# BTI Filling Tests - With extra columns
let
    times = get_bti_test_times()
    energies = get_bti_test_energies()
    
    mock_headers = Dict{String,Any}()
    mock_metadata = FITSMetadata(
        "",  # filepath
        1,   # hdu
        "keV",  # energy_units
        Dict{String,Vector}(  # extra_columns
            "pi_channel" => [100, 200, 300, 400],
            "detector_id" => [1, 2, 1, 2]
        ),
        mock_headers,  # headers
        nothing,       # gti
        nothing        # gti_source
    )
    
    el = EventList(times, energies, mock_metadata)
    gtis = [1.0 2.5; 5.5 7.5]
    
    original_length = length(el.times)
    fill_bad_time_intervals!(el, gtis, dt=1.0, random_fill_threshold=5.0,
                           rng=Random.MersenneTwister(123))
    
    @test length(el.meta.extra_columns["pi_channel"]) == length(el.times)
    @test length(el.meta.extra_columns["detector_id"]) == length(el.times)
    
    synthetic_pi = el.meta.extra_columns["pi_channel"][(original_length+1):end]
    @test all(pi -> pi ∈ [100, 200, 300, 400], synthetic_pi)
end

# GTI validation tests - Valid GTIs
let
    valid_gtis = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    @test check_gtis(valid_gtis) === true
end

# GTI validation tests - Invalid GTI dimensions
let
    @test_throws ArgumentError begin
        invalid_gtis = [1.0 2.0 3.0; 4.0 5.0 6.0]
        check_gtis(invalid_gtis)
    end
end

# GTI validation tests - Invalid time ordering
let
    @test_throws ArgumentError begin
        invalid_gtis = [2.0 1.0; 4.0 3.0]
        check_gtis(invalid_gtis)
    end
end

# BTI calculation tests - Basic BTI calculation
let
    gtis = [2.0 4.0; 6.0 8.0]
    btis = get_btis(gtis, 1.0, 9.0)
    
    expected_btis = [1.0 2.0; 4.0 6.0; 8.0 9.0]
    @test size(btis) == size(expected_btis)
    @test btis ≈ expected_btis
end

# BTI calculation tests - No BTIs (continuous GTI)
let
    gtis = [1.0 9.0]
    btis = get_btis(gtis, 1.0, 9.0)
    
    @test size(btis, 1) == 0
end

# BTI calculation tests - Empty GTIs
let
    gtis = reshape(Float64[], 0, 2)
    btis = get_btis(gtis, 1.0, 9.0)
    
    @test size(btis) == (1, 2)
    @test btis ≈ [1.0 9.0]
end
# test_join_equal_gti_boundaries
let
    gti = [1.16703354e8 1.16703386e8;
           1.16703386e8 1.16703418e8;
           1.16703418e8 1.16703450e8;
           1.16703450e8 1.16703482e8;
           1.16703482e8 1.16703514e8]
    
    newg = join_equal_gti_boundaries(gti)
    @test newg ≈ reshape([1.16703354e8 1.16703514e8], 1, 2)
end

# test_split_gtis_by_exposure
let
    gtis = [0.0 30.0; 86450.0 86460.0]
    new_gtis = split_gtis_by_exposure(gtis, 400.0)
    @test new_gtis[1] ≈ [0.0 30.0]
    @test new_gtis[2] ≈ [86450.0 86460.0]
end

# test_split_gtis_at_indices
let
    gtis = [0.0 30.0; 50.0 60.0; 80.0 90.0]
    new_gtis = split_gtis_at_indices(gtis, 1)
    @test new_gtis[1] ≈ [0.0 30.0]
    @test new_gtis[2] ≈ [50.0 60.0; 80.0 90.0]
end

# test_check_separate
let
    # Overlapping case
    gti1 = [1.0 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 13.2]
    gti2 = [2.0 5.0; 6.0 9.0; 11.4 14.0]
    @test !check_separate(gti1, gti2)  # overlapping case

    # Non-overlapping case
    gti3 = [1.0 2.0; 4.0 5.0]
    gti4 = [6.0 7.0; 8.0 9.0]
    @test check_separate(gti3, gti4)   # non-overlapping case

    # Empty case
    gti5 = Matrix{Float64}(undef, 0, 2)
    @test check_separate(gti1, gti5)
    @test check_separate(gti5, gti1)
    
    # Touching intervals (should be separate)
    gti6 = [1.0 2.0]
    gti7 = [2.0 3.0]
    @test check_separate(gti6, gti7)
    
    # Single overlapping interval
    gti8 = [1.0 3.0]
    gti9 = [2.0 4.0]
    @test !check_separate(gti8, gti9)
end

# test_append_gtis
let
    # Basic append non-overlapping GTIs
    gti1 = [1.0 2.0; 4.0 5.0]
    gti2 = [6.0 7.0; 8.0 9.0]
    @test append_gtis(gti1, gti2) ≈ [1.0 2.0; 4.0 5.0; 6.0 7.0; 8.0 9.0]
    
    # Test overlapping case - should merge, not throw error
    gti3 = [3.0 4.5; 8.0 9.0]
    expected = [1.0 2.0; 3.0 5.0; 8.0 9.0]  # [4.0,5.0] and [3.0,4.5] merge to [3.0,5.0]
    @test append_gtis(gti1, gti3) ≈ expected
    
    # Append touching GTIs
    gti4 = [1.0 2.0]
    gti5 = [2.0 3.0]
    result = append_gtis(gti4, gti5)
    @test result ≈ [1.0 3.0]
    
    # Append overlapping GTIs
    gti6 = [1.0 3.0]
    gti7 = [2.0 4.0]
    result = append_gtis(gti6, gti7)
    @test result ≈ [1.0 4.0]
end

# test_gti_border_bins
let
    # Test with simple time array
    times = collect(range(0.5, 2.5, length=3))  # [0.5, 1.5, 2.5]
    start_bins, stop_bins = gti_border_bins([0.0 2.0], times)
    @test start_bins == [0]  # First bin starts at index 0 (Python-style)
    @test stop_bins == [2]   # Last bin is at index 2

    # Test with many bins - fixed bin counting
    dt = 0.0001
    times_dense = collect(range(0.0, 2.0-dt, step=dt)) .+ dt/2  # Note the -dt to get exact number of bins
    @test length(times_dense) == 20000  # Verify the length first
    start_bins, stop_bins = gti_border_bins([0.0 2.0], times_dense)
    @test start_bins == [0]
    @test stop_bins == [20000]

    # Additional test case for boundary conditions
    times_edge = [0.5, 1.0, 1.5, 2.0]
    start_bins, stop_bins = gti_border_bins([0.75 1.75], times_edge)
    @test start_bins == [1]
    @test stop_bins == [3]
end

# test_check_gtis
let
    # Valid GTI
    gti = [1.0 2.0; 3.0 4.0]
    @test check_gtis(gti) == true
    
    # Empty GTI
    gti_empty = Matrix{Float64}(undef, 0, 2)
    @test check_gtis(gti_empty) == true
    
    # Invalid GTI - wrong number of columns
    gti_wrong_cols = [1.0 2.0 3.0]
    @test_throws ArgumentError check_gtis(gti_wrong_cols)
    
    # Invalid GTI - start >= stop
    gti_invalid = [2.0 1.0]
    @test_throws ArgumentError check_gtis(gti_invalid)
    
    # Invalid GTI - start == stop
    gti_equal = [1.0 1.0]
    @test_throws ArgumentError check_gtis(gti_equal)
end

# test_merge_overlapping_gtis
let
    # Basic merging
    gtis = [0.0 2.0; 1.0 3.0; 4.0 5.0]
    result = merge_overlapping_gtis(gtis)
    @test result ≈ [0.0 3.0; 4.0 5.0]
    
    # Touching intervals
    gtis = [0.0 1.0; 1.0 2.0; 2.0 3.0]
    result = merge_overlapping_gtis(gtis)
    @test result ≈ [0.0 3.0]
    
    # No overlap
    gtis = [0.0 1.0; 2.0 3.0; 4.0 5.0]
    result = merge_overlapping_gtis(gtis)
    @test result ≈ gtis
    
    # Single interval
    gtis = [0.0 1.0]
    result = merge_overlapping_gtis(gtis)
    @test result ≈ gtis
    
    # Empty input
    gtis = Matrix{Float64}(undef, 0, 2)
    result = merge_overlapping_gtis(gtis)
    @test isempty(result)
    
    # Complex overlapping pattern
    gtis = [1.0 3.0; 2.0 5.0; 7.0 9.0; 8.0 10.0]
    result = merge_overlapping_gtis(gtis)
    @test result ≈ [1.0 5.0; 7.0 10.0]
end

# test_join_gtis
let
    # Non-overlapping GTIs
    gti0 = [0.0 1.0; 2.0 3.0]
    gti1 = [10.0 11.0; 12.0 13.0]
    result = join_gtis(gti0, gti1)
    expected = [0.0 1.0; 2.0 3.0; 10.0 11.0; 12.0 13.0]
    @test result ≈ expected

    # Overlapping GTIs
    gti0 = [0.0 1.0; 2.0 3.0; 4.0 8.0]
    gti1 = [7.0 8.0; 10.0 11.0; 12.0 13.0]
    result = join_gtis(gti0, gti1)
    expected = [0.0 1.0; 2.0 3.0; 4.0 8.0; 10.0 11.0; 12.0 13.0]
    @test result ≈ expected

    # GTIs meeting in middle
    gti0 = [0.0 1.0; 2.0 3.0; 4.0 8.0]
    gti1 = [1.0 2.0; 3.0 4.0]
    result = join_gtis(gti0, gti1)
    expected = [0.0 8.0]
    @test result ≈ expected
    
    # Empty GTI cases
    gti_empty = Matrix{Float64}(undef, 0, 2)
    gti_data = [1.0 2.0; 3.0 4.0]
    @test join_gtis(gti_empty, gti_data) ≈ gti_data
    @test join_gtis(gti_data, gti_empty) ≈ gti_data
    @test isempty(join_gtis(gti_empty, gti_empty))
end

# test_cross_two_gtis
let
    # Simple intersection
    gti1 = [1.0 4.0]
    gti2 = [2.0 3.0]
    result = cross_two_gtis(gti1, gti2)
    @test result ≈ [2.0 3.0]
    
    # Multiple intervals
    gti1 = [1.0 5.0; 7.0 10.0]
    gti2 = [2.0 4.0; 8.0 9.0]
    result = cross_two_gtis(gti1, gti2)
    @test result ≈ [2.0 4.0; 8.0 9.0]
    
    # No intersection
    gti1 = [1.0 2.0]
    gti2 = [3.0 4.0]
    result = cross_two_gtis(gti1, gti2)
    @test isempty(result)
    
    # Empty GTI
    gti_empty = Matrix{Float64}(undef, 0, 2)
    gti_data = [1.0 2.0]
    @test isempty(cross_two_gtis(gti_empty, gti_data))
    @test isempty(cross_two_gtis(gti_data, gti_empty))
    
    # Partial overlap
    gti1 = [1.0 3.0]
    gti2 = [2.0 4.0]
    result = cross_two_gtis(gti1, gti2)
    @test result ≈ [2.0 3.0]
end

# test_cross_gtis
let
    # Basic intersection
    gti1 = [1.0 4.0]
    gti2 = [2.0 5.0]
    result = cross_gtis([gti1, gti2])
    @test result ≈ [2.0 4.0]

    # Complex intersection
    gti1 = [1.0 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 13.2]
    gti2 = [2.0 5.0; 6.0 9.0; 11.4 14.0]
    result = cross_gtis([gti1, gti2])
    @test result ≈ [4.0 5.0; 7.0 9.0; 12.2 13.2]

    # Single GTI
    gti1 = [1.0 2.0; 4.0 5.0; 7.0 10.0]
    result = cross_gtis([gti1])
    @test result ≈ gti1

    # Empty intersection
    gti1 = [2.0 3.0]
    gti2 = [3.0 4.0]
    result = cross_gtis([gti1, gti2])
    @test isempty(result)
    
    gti3 = [3.0 5.0]
    result = cross_gtis([gti1, gti2, gti3])
    @test isempty(result)

    # Wide GTI intersection
    gti1 = [1.0 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 13.2]
    gti2 = [0.5 14.0]
    result1 = cross_gtis([gti1, gti2])
    result2 = cross_gtis([gti2, gti1])
    @test result1 ≈ gti1
    @test result2 ≈ gti1

    # Partial overlap
    gti1 = [1.5 12.5]
    gti2 = [1.0 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 13.2]
    result1 = cross_gtis([gti1, gti2])
    result2 = cross_gtis([gti2, gti1])
    expected = [1.5 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 12.5]
    @test result1 ≈ expected
    @test result2 ≈ expected

    # Complex overlap
    gti1 = [1.0 2.0; 4.0 5.0; 7.0 10.0; 11.0 11.2; 12.2 13.2]
    gti2 = [0.5 3.0; 4.5 4.7; 10.0 14.0]
    result1 = cross_gtis([gti1, gti2])
    result2 = cross_gtis([gti2, gti1])
    expected = [1.0 2.0; 4.5 4.7; 11.0 11.2; 12.2 13.2]
    @test result1 ≈ expected
    @test result2 ≈ expected
end

# test_merge_gtis
let
    # Setup test GTIs
    gti1 = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    gti2 = [1.0 2.0]
    gti3 = [2.0 3.0]
    gti4 = [4.0 5.0]

    # Test empty GTIs for all methods
    for method in ["intersection", "union", "infer", "append", "none"]
        @test isnothing(merge_gtis(Matrix{Float64}[], method))
        @test isnothing(merge_gtis([Matrix{Float64}(undef,0,2)], method))
        @test isnothing(merge_gtis([Matrix{Float64}(undef,0,2), Matrix{Float64}(undef,0,2)], method))
    end

    # Test single GTI for all methods
    for method in ["intersection", "union", "infer", "append"]
        @test merge_gtis([gti1], method) ≈ gti1
    end

    # Test "none" method
    @test merge_gtis([gti1], "none") ≈ [1.0 6.0]
    @test merge_gtis([gti1, gti2], "none") ≈ [1.0 6.0]

    # Test intersection method
    @test merge_gtis([gti1, gti2], "intersection") ≈ [1.0 2.0]
    @test isnothing(merge_gtis([gti1, gti2, gti3], "intersection"))
    @test isnothing(merge_gtis([gti2, gti3], "intersection"))

    # Test union method
    @test merge_gtis([gti1, gti2], "union") ≈ gti1
    @test merge_gtis([gti1, gti3], "union") ≈ [1.0 4.0; 5.0 6.0]

    # Test append method
    @test merge_gtis([gti2, gti4], "append") ≈ [1.0 2.0; 4.0 5.0]
    @test merge_gtis([gti2, gti3], "append") ≈ [1.0 3.0]

    # Test infer method
    @test merge_gtis([gti1, gti2], "infer") ≈ [1.0 2.0]
    @test merge_gtis([gti2, gti4], "infer") ≈ [1.0 2.0; 4.0 5.0]
end