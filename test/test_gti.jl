# Helper function to create mock EventList for testing
function create_test_eventlist(times::Vector{Float64}, energies::Union{Vector{Float64}, Nothing}=nothing)
    mock_headers = Dict{String,Any}()
    mock_metadata = FITSMetadata(
        "", 1, "keV",
        Dict{String,Vector}(),
        mock_headers
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
    @test_throws ArgumentError Stingray.check_gtis([[0 1 4]; [0 3 4];]) 
end

# test_check_gtis_values
let
    @test_throws ArgumentError Stingray.check_gtis([[0 2]; [1 3];])

    @test_throws ArgumentError Stingray.check_gtis([[1 0];])
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
        "", 1, "keV",
        Dict{String,Vector}(
            "pi_channel" => [100, 200, 300, 400],
            "detector_id" => [1, 2, 1, 2]
        ),
        mock_headers
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
    @test check_gtis(valid_gtis) === nothing
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

# Error calculation tests
let
    counts = [0, 1, 4, 9, 16]
    
    poisson_errors = calculate_errors(counts, :poisson)
    @test poisson_errors ≈ [1.0, 1.0, 2.0, 3.0, 4.0]
    
    provided_errors = [0.5, 1.2, 2.1, 3.3, 4.0]
    gaussian_errors = calculate_errors(counts, :gaussian, gaussian_errors=provided_errors)
    @test gaussian_errors ≈ provided_errors
end

# Error calculation tests - Invalid method
let
    @test_throws ArgumentError begin
        counts = [1, 2, 3]
        calculate_errors(counts, :invalid_method)
    end
end

# Error calculation tests - Gaussian without errors
let
    @test_throws ArgumentError begin
        counts = [1, 2, 3]
        calculate_errors(counts, :gaussian)
    end
end

# Error calculation tests - Gaussian length mismatch
let
    @test_throws ArgumentError begin
        counts = [1, 2, 3]
        wrong_length_errors = [1.0, 2.0]
        calculate_errors(counts, :gaussian, gaussian_errors=wrong_length_errors)
    end
end

# EventList energy check tests
let
    times = [1.0, 2.0, 3.0]
    energies = [10.0, 20.0, 30.0]
    el_with_energies = create_test_eventlist(times, energies)
    @test has_energies(el_with_energies) == true
    
    el_without_energies = create_test_eventlist(times, nothing)
    @test has_energies(el_without_energies) == false
end