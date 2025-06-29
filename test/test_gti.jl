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