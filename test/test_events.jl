

# Test data path (if using test data directory)
# TEST_DATA_PATH = joinpath(@__DIR__, "_data", "testdata")

# Helper function to generate mock data
function mock_data(times, energies; energy_column = "ENERGY")
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "sample.fits")

    # Create a sample FITS file
    FITS(sample_file, "w") do f
        # Create primary HDU with a small array instead of empty
        # Use a single element array instead of empty
        write(f, [0])

        # Create event table in HDU 2
        table = Dict{String,Array}()
        table["TIME"] = times
        table[energy_column] = energies
        table["INDEX"] = collect(1:length(times))
        write(f, table)
    end
    sample_file
end

# Test basic EventList creation and validation
let
    # Test valid construction with simplified constructor
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 15.0, 25.0, 30.0]

    ev = EventList(times, energies)
    @test ev.times == times
    @test ev.energies == energies
    @test length(ev) == 5
    @test has_energies(ev)

    # Test with no energies
    ev_no_energy = EventList(times)
    @test ev_no_energy.times == times
    @test isnothing(ev_no_energy.energies)
    @test !has_energies(ev_no_energy)

end

# Test accessor functions
let
    times_vec = [1.0, 2.0, 3.0]
    energies_vec = [10.0, 20.0, 30.0]

    ev = EventList(times_vec, energies_vec)

    # Test times() accessor
    @test times(ev) === ev.times
    @test times(ev) == times_vec

    # Test energies() accessor  
    @test energies(ev) === ev.energies
    @test energies(ev) == energies_vec

    # Test with nothing energies
    ev_no_energy = EventList(times_vec)
    @test isnothing(energies(ev_no_energy))

end

# Test Base interface methods
let
    times_vec = [1.0, 2.0, 3.0, 4.0]

    ev = EventList(times_vec)

    # Test length
    @test length(ev) == 4
    @test length(ev) == length(times_vec)

    # Test size
    @test size(ev) == (4,)
    @test size(ev) == (length(times_vec),)

end

# Test filter_time! function (in-place filtering by time)
let
    # Test basic time filtering - keep times >= 3.0
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]

    ev = EventList(times, energies)
    filter_time!(t -> t >= 3.0, ev)

    @test ev.times == [3.0, 4.0, 5.0]
    @test ev.energies == [30.0, 40.0, 50.0]
    @test length(ev) == 3

    # Test filtering that removes all elements
    ev_empty = EventList([1.0, 2.0], [10.0, 20.0])
    filter_time!(t -> t > 10.0, ev_empty)
    @test length(ev_empty) == 0
    @test ev_empty.times == Float64[]
    @test ev_empty.energies == Float64[]

    # Test filtering with no energies
    ev_no_energy = EventList([1.0, 2.0, 3.0, 4.0])
    filter_time!(t -> t > 2.0, ev_no_energy)
    @test ev_no_energy.times == [3.0, 4.0]
    @test isnothing(ev_no_energy.energies)
    @test length(ev_no_energy) == 2

    # Test filtering with extra columns
    times_extra = [1.0, 2.0, 3.0, 4.0]
    energies_extra = [10.0, 20.0, 30.0, 40.0]
    dummy_meta = FITSMetadata{Dict{String,Any}}(
        "",
        1,
        nothing,
        Dict("INDEX" => [1, 2, 3, 4]),
        Dict{String,Any}(),
    )
    ev_extra = EventList(times_extra, energies_extra, dummy_meta)

    filter_time!(t -> t >= 2.5, ev_extra)
    @test ev_extra.times == [3.0, 4.0]
    @test ev_extra.energies == [30.0, 40.0]
    @test ev_extra.meta.extra_columns["INDEX"] == [3, 4]

end

# Test filter_energy! function (in-place filtering by energy)
let
    # Test basic energy filtering - keep energies >= 25.0
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]

    ev = EventList(times, energies)
    filter_energy!(e -> e >= 25.0, ev)

    @test ev.times == [3.0, 4.0, 5.0]
    @test ev.energies == [30.0, 40.0, 50.0]
    @test length(ev) == 3

    # Test filtering that removes all elements
    ev_all_removed = EventList([1.0, 2.0], [10.0, 20.0])
    filter_energy!(e -> e > 100.0, ev_all_removed)
    @test length(ev_all_removed) == 0
    @test ev_all_removed.times == Float64[]
    @test ev_all_removed.energies == Float64[]

    # Test error when no energies present
    ev_no_energy = EventList([1.0, 2.0, 3.0])
    @test_throws AssertionError filter_energy!(e -> e > 10.0, ev_no_energy)

    # Test filtering with extra columns
    times_extra = [1.0, 2.0, 3.0, 4.0]
    energies_extra = [15.0, 25.0, 35.0, 45.0]
    dummy_meta = FITSMetadata{Dict{String,Any}}(
        "",
        1,
        nothing,
        Dict("DETX" => [0.1, 0.2, 0.3, 0.4]),
        Dict{String,Any}(),
    )
    ev_extra = EventList(times_extra, energies_extra, dummy_meta)

    filter_energy!(e -> e >= 30.0, ev_extra)
    @test ev_extra.times == [3.0, 4.0]
    @test ev_extra.energies == [35.0, 45.0]
    @test ev_extra.meta.extra_columns["DETX"] == [0.3, 0.4]

    println("✓ filter_energy! function tests passed")
end

# Test filter_on! function (generic in-place filtering)
let
    # Test filtering on times using filter_on!
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]

    ev = EventList(times, energies)
    filter_on!(t -> t % 2 == 0, ev.times, ev)  # Keep even time values (2.0, 4.0)

    @test ev.times == [2.0, 4.0]
    @test ev.energies == [20.0, 40.0]
    @test length(ev) == 2

    # Test filtering on energies using filter_on!
    times2 = [1.0, 2.0, 3.0, 4.0]
    energies2 = [15.0, 25.0, 35.0, 45.0]
    ev2 = EventList(times2, energies2)
    filter_on!(e -> e > 30.0, ev2.energies, ev2)  # Keep energies > 30

    @test ev2.times == [3.0, 4.0]
    @test ev2.energies == [35.0, 45.0]

    # Test assertion error for mismatched sizes
    times3 = [1.0, 2.0, 3.0]
    energies3 = [10.0, 20.0, 30.0]
    ev3 = EventList(times3, energies3)
    wrong_size_col = [1.0, 2.0]  # Wrong size
    @test_throws AssertionError filter_on!(x -> x > 1.0, wrong_size_col, ev3)

    # Test with extra columns
    times_extra = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies_extra = [10.0, 20.0, 30.0, 40.0, 50.0]
    dummy_meta = FITSMetadata{Dict{String,Any}}(
        "",
        1,
        nothing,
        Dict("FLAG" => [1, 0, 1, 0, 1]),
        Dict{String,Any}(),
    )
    ev_extra = EventList(times_extra, energies_extra, dummy_meta)

    # Filter based on FLAG column (keep where FLAG == 1)
    filter_on!(flag -> flag == 1, ev_extra.meta.extra_columns["FLAG"], ev_extra)
    @test ev_extra.times == [1.0, 3.0, 5.0]
    @test ev_extra.energies == [10.0, 30.0, 50.0]
    @test ev_extra.meta.extra_columns["FLAG"] == [1, 1, 1]

    println("✓ filter_on! function tests passed")
end

# Test non-mutating filter functions (filter_time and filter_energy)
let
    # Test filter_time (non-mutating)
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]

    ev_original = EventList(times, energies)
    ev_filtered = filter_time(t -> t >= 3.0, ev_original)

    # Original should be unchanged
    @test ev_original.times == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test ev_original.energies == [10.0, 20.0, 30.0, 40.0, 50.0]
    @test length(ev_original) == 5

    # Filtered should have new values
    @test ev_filtered.times == [3.0, 4.0, 5.0]
    @test ev_filtered.energies == [30.0, 40.0, 50.0]
    @test length(ev_filtered) == 3

    # Test filter_energy (non-mutating)
    ev_original2 = EventList(times, energies)
    ev_filtered2 = filter_energy(e -> e <= 30.0, ev_original2)

    # Original should be unchanged
    @test ev_original2.times == times
    @test ev_original2.energies == energies

    # Filtered should have new values
    @test ev_filtered2.times == [1.0, 2.0, 3.0]
    @test ev_filtered2.energies == [10.0, 20.0, 30.0]
    @test length(ev_filtered2) == 3

    # Test with no energies
    ev_no_energy = EventList([1.0, 2.0, 3.0, 4.0])
    ev_filtered_no_energy = filter_time(t -> t > 2.0, ev_no_energy)

    @test ev_no_energy.times == [1.0, 2.0, 3.0, 4.0]  # Original unchanged
    @test ev_filtered_no_energy.times == [3.0, 4.0]
    @test isnothing(ev_filtered_no_energy.energies)

    # Test filter_energy error with no energies
    @test_throws AssertionError filter_energy(e -> e > 10.0, ev_no_energy)

    println("✓ Non-mutating filter function tests passed")
end

# Test complex filtering scenarios
let
    # Test multiple sequential filters
    times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    energies = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

    ev = EventList(times, energies)

    # First filter by time (keep t >= 3.0)
    filter_time!(t -> t >= 3.0, ev)
    @test length(ev) == 6
    @test ev.times == [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    @test ev.energies == [20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

    # Then filter by energy (keep e <= 35.0)
    filter_energy!(e -> e <= 35.0, ev)
    @test length(ev) == 4
    @test ev.times == [3.0, 4.0, 5.0, 6.0]
    @test ev.energies == [20.0, 25.0, 30.0, 35.0]

    # Test edge case: empty result after filtering
    ev_edge = EventList([1.0, 2.0], [10.0, 20.0])
    filter_time!(t -> t > 5.0, ev_edge)
    @test length(ev_edge) == 0
    @test ev_edge.times == Float64[]
    @test ev_edge.energies == Float64[]

    # Test edge case: no filtering needed (all pass)
    ev_all_pass = EventList([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])
    original_times = copy(ev_all_pass.times)
    original_energies = copy(ev_all_pass.energies)
    filter_time!(t -> t > 0.0, ev_all_pass)  # All should pass
    @test ev_all_pass.times == original_times
    @test ev_all_pass.energies == original_energies
    @test length(ev_all_pass) == 3

    println("✓ Complex filtering scenario tests passed")
end

# Test readevents basic functionality with mock data
let
    mock_times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
    mock_energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
    sample = mock_data(mock_times, mock_energies)

    # Try reading mock data
    data = readevents(sample)
    @test data.times == mock_times
    @test data.energies == mock_energies
    @test length(data.meta.headers) >= 1  # Should have at least primary header
    @test length(data.meta.extra_columns) == 0

    # Test reading with extra columns
    data = readevents(sample; extra_columns = ["INDEX"])
    @test length(data.meta.extra_columns) == 1
    @test haskey(data.meta.extra_columns, "INDEX")
    @test data.meta.extra_columns["INDEX"] == collect(1:5)

    println("✓ Basic readevents tests passed")
end

# Test readevents HDU handling
let
    # Test with events in HDU 3 instead of default HDU 2
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "hdu3_sample.fits")
    f = FITS(sample_file, "w")
    write(f, [0])  # Primary HDU with non-empty array

    # Empty table in HDU 2
    empty_table = Dict{String,Array}()
    empty_table["OTHER"] = Float64[1.0, 2.0]
    write(f, empty_table)

    # Event data in HDU 3
    times = Float64[1.0, 2.0, 3.0]
    energies = Float64[10.0, 20.0, 30.0]
    event_table = Dict{String,Array}()
    event_table["TIME"] = times
    event_table["ENERGY"] = energies
    write(f, event_table)
    close(f)

    # Test specifying specific HDU
    data_hdu3 = readevents(sample_file, hdu = 3)
    @test data_hdu3.times == times
    @test data_hdu3.energies == energies

end

# Test readevents alternative energy columns
let
    # Test with PI column
    times = Float64[1.0, 2.0, 3.0]
    pi_values = Float64[100.0, 200.0, 300.0]

    pi_file = mock_data(times, pi_values; energy_column = "PI")

    data = readevents(pi_file)
    @test data.times == times
    @test data.energies == pi_values
    @test data.meta.energy_units == "PI"

    # Test with PHA column
    pha_file = mock_data(times, pi_values; energy_column = "PHA")

    data_pha = readevents(pha_file)
    @test data_pha.times == times
    @test data_pha.energies == pi_values
    @test data_pha.meta.energy_units == "PHA"

end

# Test readevents missing columns
let
    # File with only TIME column
    test_dir = mktempdir()
    time_only_file = joinpath(test_dir, "time_only.fits")
    f = FITS(time_only_file, "w")
    write(f, [0])  # Non-empty primary HDU

    times = Float64[1.0, 2.0, 3.0]
    table = Dict{String,Array}()
    table["TIME"] = times
    write(f, table)
    close(f)

    data = readevents(time_only_file)
    @test data.times == times
    @test isnothing(data.energies)
    @test isnothing(data.meta.energy_units)

end

# Test error handling
let
    # Test non-existent file
    @test_throws Exception readevents("non_existent_file.fits")

    # Test invalid FITS file
    test_dir = mktempdir()
    invalid_file = joinpath(test_dir, "invalid.fits")
    open(invalid_file, "w") do io
        write(io, "This is not a FITS file")
    end
    @test_throws Exception readevents(invalid_file)

end

# Test case insensitive column names
let
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "case_test.fits")

    f = FITS(sample_file, "w")

    # Create primary HDU with valid data
    primary_data = reshape([1.0], 1, 1)
    write(f, primary_data)

    # Use lowercase column names
    times = Float64[1.0, 2.0, 3.0]
    energies = Float64[10.0, 20.0, 30.0]

    table = Dict{String,Array}()
    table["time"] = times  # lowercase
    table["energy"] = energies  # lowercase
    write(f, table)
    close(f)

    data = readevents(sample_file)
    @test data.times == times
    @test data.energies == energies

end
