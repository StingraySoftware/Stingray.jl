"""
    create_synthetic_events(n_events::Int=1000; mission::String="nustar", 
                           time_range::Tuple{Float64,Float64}=(0.0, 1000.0),
                           pi_range::Tuple{Int,Int}=(1, 4096),
                           T::Type=Float64) -> EventList{Vector{T}, FITSMetadata}

Create synthetic X-ray event data for testing purposes.

This function generates synthetic X-ray event data with proper time ordering,
energy calibration, and mission-specific metadata. The generated events include:
- Random event times within the specified range
- Random PI (Pulse Invariant) channels
- Energy values calculated using mission-specific calibration
- Detector ID assignments
- Realistic FITS-style metadata headers

# Arguments
- `n_events::Int=1000`: Number of events to generate
- `mission::String="nustar"`: Mission name for calibration and metadata
- `time_range::Tuple{Float64,Float64}=(0.0, 1000.0)`: Time range for events (start, stop)
- `pi_range::Tuple{Int,Int}=(1, 4096)`: PI channel range
- `T::Type=Float64`: Numeric type for the data

# Returns
- `EventList{Vector{T}, FITSMetadata}`: Synthetic event list with times, energies, and metadata

# Examples
```julia
# Default NuSTAR events
events = create_synthetic_events()

# Custom XMM events
events = create_synthetic_events(5000; mission="xmm", time_range=(100.0, 200.0))

# NICER events with custom PI range
events = create_synthetic_events(2000; mission="nicer", pi_range=(50, 1000))
```
"""
function create_synthetic_events(n_events::Int=1000; 
                               mission::String="nustar",
                               time_range::Tuple{Float64,Float64}=(0.0, 1000.0),
                               pi_range::Tuple{Int,Int}=(1, 4096),
                               T::Type=Float64)
    
    # Input validation
    if n_events <= 0
        throw(ArgumentError("Number of events must be positive"))
    end
    
    if time_range[1] >= time_range[2]
        throw(ArgumentError("Time range start must be less than end"))
    end
    
    if pi_range[1] > pi_range[2]
        throw(ArgumentError("PI range start must be less than or equal to end"))
    end
    
    # Get mission support for calibration
    mission_support = get_mission_support(mission)
    
    # Generate random event times and sort them
    times = sort(rand(n_events) * (time_range[2] - time_range[1]) .+ time_range[1])
    times = convert(Vector{T}, times)
    
    # Generate random PI channels
    pi_channels = rand(pi_range[1]:pi_range[2], n_events)
    
    # Apply mission-specific calibration to get energies
    energies = convert(Vector{T}, apply_calibration(mission_support, pi_channels))
    
    # Generate detector IDs (assume 4-detector system like NuSTAR)
    detector_ids = rand(0:3, n_events)
    
    # Create extra columns dictionary
    extra_columns = Dict{String, Vector}(
        "PI" => pi_channels,
        "DETID" => detector_ids
    )
    
    # Create synthetic filename
    filepath = "synthetic_$(mission).fits"
    
    # Create FITS headers using FITSIO.FITSHeader structure
    # Note: This is a simplified header for testing - in real usage, 
    # headers would come from actual FITS files
    header_dict = Dict{String, Any}(
        "TELESCOP" => uppercase(mission),
        "INSTRUME" => "TEST_INSTRUMENT",
        "OBJECT" => "TEST_SOURCE",
        "RA_NOM" => 180.0,
        "DEC_NOM" => 0.0,
        "EQUINOX" => 2000.0,
        "RADECSYS" => "FK5",
        "TSTART" => time_range[1],
        "TSTOP" => time_range[2],
        "EXPOSURE" => time_range[2] - time_range[1],
        "ONTIME" => time_range[2] - time_range[1],
        "LIVETIME" => time_range[2] - time_range[1],
        "NAXIS2" => n_events,
        "TIMESYS" => "TT",
        "TIMEREF" => "LOCAL",
        "TIMEUNIT" => "s",
        "MJDREFI" => 55197,  # Standard MJD reference for many missions
        "MJDREFF" => 0.00076601852,
        "CLOCKCORR" => "T",
        "DATE-OBS" => "2020-01-01T00:00:00",
        "DATE-END" => "2020-01-01T01:00:00"
    )
    
    # Add mission-specific header information
    if lowercase(mission) == "nustar"
        header_dict["DETNAM"] = "TEST_DET"
        header_dict["PIFLTCOR"] = "T"
    elseif lowercase(mission) == "xmm"
        header_dict["FILTER"] = "Medium"
        header_dict["SUBMODE"] = "PrimeFullWindow"
    elseif lowercase(mission) == "nicer"
        header_dict["DETNAM"] = "TEST_NICER"
        header_dict["FILTFILE"] = "NONE"
    elseif lowercase(mission) in ["chandra", "axaf"]
        header_dict["DETNAM"] = "ACIS-S"
        header_dict["GRATING"] = "NONE"
    elseif lowercase(mission) == "xte"
        header_dict["DETNAM"] = "PCA"
        header_dict["LAYERS"] = "ALL"
    elseif lowercase(mission) == "ixpe"
        header_dict["DETNAM"] = "GPD"
        header_dict["POLMODE"] = "ON"
    end
    
    # Create a simple header structure (for testing purposes)
    # In real usage, this would be a proper FITSIO.FITSHeader
    headers = header_dict
    
    # Create FITSMetadata
    metadata = FITSMetadata(
        filepath,           # filepath
        2,                  # hdu (typical for event data)
        "ENERGY",          # energy_units (after calibration)
        extra_columns,     # extra_columns
        headers            # headers
    )
    
    # Create and return EventList
    return EventList(times, energies, metadata)
end
# Test: Basic Synthetic Event Creation - Default Parameters
let
    events = create_synthetic_events()
    
    # Test basic structure
    @test isa(events, EventList)
    @test length(events) == 1000  # Default n_events
    @test length(times(events)) == 1000
    @test length(energies(events)) == 1000
    @test events.meta.filepath == "synthetic_nustar.fits"
    
    # Test times are sorted and in range
    @test issorted(times(events))
    @test all(0.0 .<= times(events) .<= 1000.0)
    
    # Test energies are positive (after calibration)
    @test all(energies(events) .> 0)
    
    # Test extra columns
    @test haskey(events.meta.extra_columns, "PI")
    @test haskey(events.meta.extra_columns, "DETID")
    @test length(events.meta.extra_columns["PI"]) == 1000
    @test length(events.meta.extra_columns["DETID"]) == 1000
    
    # Test PI channels are in expected range
    @test all(1 .<= events.meta.extra_columns["PI"] .<= 4096)
    
    # Test detector IDs are in expected range
    @test all(0 .<= events.meta.extra_columns["DETID"] .<= 3)
    
    # Test metadata structure
    @test isa(events.meta, FITSMetadata)
    @test events.meta.hdu == 2
    @test events.meta.energy_units == "ENERGY"
    
    # Test metadata content
    @test events.meta.headers["TELESCOP"] == "NUSTAR"
    @test events.meta.headers["INSTRUME"] == "TEST_INSTRUMENT"
    @test events.meta.headers["NAXIS2"] == 1000
    @test events.meta.headers["TSTART"] == 0.0
    @test events.meta.headers["TSTOP"] == 1000.0
end

# Test: Basic Synthetic Event Creation - Custom Parameters
let
    n_events = 500
    mission = "xmm"
    time_range = (100.0, 200.0)
    pi_range = (50, 1000)
    
    events = create_synthetic_events(n_events; 
                                   mission=mission, 
                                   time_range=time_range, 
                                   pi_range=pi_range)
    
    # Test custom parameters are respected
    @test length(events) == n_events
    @test length(times(events)) == n_events
    @test length(energies(events)) == n_events
    @test events.meta.filepath == "synthetic_xmm.fits"
    
    # Test time range
    @test all(time_range[1] .<= times(events) .<= time_range[2])
    @test events.meta.headers["TSTART"] == time_range[1]
    @test events.meta.headers["TSTOP"] == time_range[2]
    
    # Test PI range
    @test all(pi_range[1] .<= events.meta.extra_columns["PI"] .<= pi_range[2])
    
    # Test mission-specific metadata
    @test events.meta.headers["TELESCOP"] == "XMM"
end

# Test: Different Missions
let
    missions = ["nustar", "xmm", "nicer", "ixpe", "axaf", "chandra", "xte"]
    
    for mission in missions
        events = create_synthetic_events(100; mission=mission)
        
        @test events.meta.filepath == "synthetic_$(mission).fits"
        @test events.meta.headers["TELESCOP"] == uppercase(mission)
        @test length(events) == 100
        @test length(times(events)) == 100
        @test length(energies(events)) == 100
        
        # Test that calibration was applied correctly
        ms = get_mission_support(mission)
        expected_energies = apply_calibration(ms, events.meta.extra_columns["PI"])
        @test energies(events) ≈ expected_energies
    end
end

# Test: Time Ordering
let
    for _ in 1:10  # Test multiple times due to randomness
        events = create_synthetic_events(100)
        @test issorted(times(events))
        
        # Test no duplicate times (very unlikely but possible)
        @test length(unique(times(events))) >= 95  # Allow some duplicates due to floating point
    end
end

# Test: Energy Calibration Consistency
let
    missions = ["nustar", "xmm", "nicer", "ixpe"]
    
    for mission in missions
        events = create_synthetic_events(200; mission=mission)
        
        # Manually verify calibration
        ms = get_mission_support(mission)
        expected_energies = apply_calibration(ms, events.meta.extra_columns["PI"])
        @test energies(events) ≈ expected_energies
        
        # Test energy ranges are reasonable for each mission
        if mission == "nustar"
            # NuSTAR: pi * 0.04 + 1.62, PI range 1-4096
            @test minimum(energies(events)) >= 1.66  # 1*0.04 + 1.62
            @test maximum(energies(events)) <= 165.46  # 4096*0.04 + 1.62
        elseif mission == "xmm"
            # XMM: pi * 0.001, PI range 1-4096
            @test minimum(energies(events)) >= 0.001
            @test maximum(energies(events)) <= 4.096
        elseif mission == "nicer"
            # NICER: pi * 0.01, PI range 1-4096
            @test minimum(energies(events)) >= 0.01
            @test maximum(energies(events)) <= 40.96
        elseif mission == "ixpe"
            # IXPE: pi / 375 * 15, PI range 1-4096
            @test minimum(energies(events)) >= 15.0/375  # ≈ 0.04
            @test maximum(energies(events)) <= 4096*15.0/375  # ≈ 163.84
        end
    end
end

# Test: Statistical Properties
let
    events = create_synthetic_events(10000)  # Large sample for statistics
    
    # Test time distribution (should be roughly uniform)
    time_hist = fit(Histogram, times(events), 0:100:1000)
    counts = time_hist.weights
    # Expect roughly equal counts in each bin (within statistical fluctuation)
    expected_count = 10000 / 10  # 10 bins
    @test all(abs.(counts .- expected_count) .< 3 * sqrt(expected_count))  # 3-sigma test
    
    # Test PI distribution (should be roughly uniform over discrete range)
    pi_hist = fit(Histogram, events.meta.extra_columns["PI"], 1:100:4096)
    pi_counts = pi_hist.weights
    # More lenient test due to discrete uniform distribution
    @test std(pi_counts) / mean(pi_counts) < 0.2  # Coefficient of variation < 20%
    
    # Test detector ID distribution
    detid_counts = [count(==(i), events.meta.extra_columns["DETID"]) for i in 0:3]
    @test all(abs.(detid_counts .- 2500) .< 3 * sqrt(2500))  # Each detector ~2500 events
end

# Test: Small Event Counts
let
    for n in [1, 2, 5, 10]
        events = create_synthetic_events(n)
        @test length(events) == n
        @test length(times(events)) == n
        @test length(energies(events)) == n
        @test length(events.meta.extra_columns["PI"]) == n
        @test length(events.meta.extra_columns["DETID"]) == n
        @test issorted(times(events))
    end
end

# Test: Large Event Counts
let
    events = create_synthetic_events(100000)
    @test length(events) == 100000
    @test length(times(events)) == 100000
    @test length(energies(events)) == 100000
    @test issorted(times(events))
    @test events.meta.headers["NAXIS2"] == 100000
end

# Test: Extreme Time Ranges
let
    # Very short time range
    events = create_synthetic_events(100; time_range=(0.0, 0.1))
    @test all(0.0 .<= times(events) .<= 0.1)
    @test events.meta.headers["TSTART"] == 0.0
    @test events.meta.headers["TSTOP"] == 0.1
    
    # Very long time range
    events = create_synthetic_events(100; time_range=(0.0, 1e6))
    @test all(0.0 .<= times(events) .<= 1e6)
    @test events.meta.headers["TSTART"] == 0.0
    @test events.meta.headers["TSTOP"] == 1e6
    
    # Negative time range
    events = create_synthetic_events(100; time_range=(-1000.0, -500.0))
    @test all(-1000.0 .<= times(events) .<= -500.0)
    @test issorted(times(events))
end

# Test: Extreme PI Ranges
let
    # Small PI range
    events = create_synthetic_events(100; pi_range=(100, 110))
    @test all(100 .<= events.meta.extra_columns["PI"] .<= 110)
    
    # Single PI value
    events = create_synthetic_events(100; pi_range=(500, 500))
    @test all(events.meta.extra_columns["PI"] .== 500)
    @test all(energies(events) .== energies(events)[1])  # All same energy
    
    # Large PI range
    events = create_synthetic_events(100; pi_range=(1, 10000))
    @test all(1 .<= events.meta.extra_columns["PI"] .<= 10000)
end

# Test: Unknown Mission
let
    # Should still work but with warning
    @test_logs (:warn, r"Mission unknown_mission not recognized") begin
        events = create_synthetic_events(100; mission="unknown_mission")
        @test events.meta.filepath == "synthetic_unknown_mission.fits"
        @test events.meta.headers["TELESCOP"] == "UNKNOWN_MISSION"
        @test length(events) == 100
        # Should use identity calibration
        @test energies(events) == Float64.(events.meta.extra_columns["PI"])
    end
end

# Test: No Missing Data
let
    events = create_synthetic_events(1000)
    
    # Check no NaN or missing values
    @test all(isfinite.(times(events)))
    @test all(isfinite.(energies(events)))
    @test all(isfinite.(events.meta.extra_columns["PI"]))
    @test all(isfinite.(events.meta.extra_columns["DETID"]))
    
    # Check no negative energies (after calibration)
    @test all(energies(events) .>= 0)
end

# Test: Correct Data Types
let
    events = create_synthetic_events(100)
    
    @test eltype(times(events)) == Float64
    @test eltype(energies(events)) == Float64
    @test eltype(events.meta.extra_columns["PI"]) <: Integer
    @test eltype(events.meta.extra_columns["DETID"]) <: Integer
    @test isa(events.meta, FITSMetadata)
    @test isa(events.meta.filepath, String)
end

# Test: Array Length Consistency
let
    for n in [10, 100, 1000, 5000]
        events = create_synthetic_events(n)
        
        @test length(events) == n
        @test length(times(events)) == n
        @test length(energies(events)) == n
        @test length(events.meta.extra_columns["PI"]) == n
        @test length(events.meta.extra_columns["DETID"]) == n
        @test events.meta.headers["NAXIS2"] == n
    end
end

# Test: Mission Name Handling
let
    # Test case sensitivity
    missions = ["NUSTAR", "nustar", "NuSTAR", "NuStar"]
    for mission in missions
        events = create_synthetic_events(100; mission=mission)
        @test events.meta.headers["TELESCOP"] == "NUSTAR"
        @test events.meta.filepath == "synthetic_$(mission).fits"
    end
end

# Test: Calibration Differences
let
    # Same PI values should give different energies for different missions
    test_pi_range = (1000, 1000)  # Fixed PI value
    
    results = Dict{String, Float64}()
    for mission in ["nustar", "xmm", "nicer", "ixpe"]
        events = create_synthetic_events(100; mission=mission, pi_range=test_pi_range)
        results[mission] = energies(events)[1]  # All should be same since PI is fixed
    end
    
    # Different missions should give different energies
    missions = collect(keys(results))
    for i in 1:length(missions)
        for j in (i+1):length(missions)
            @test results[missions[i]] ≠ results[missions[j]]
        end
    end
end

# Test: Metadata Consistency
let
    missions = ["nustar", "xmm", "nicer", "ixpe", "axaf", "chandra"]
    
    for mission in missions
        events = create_synthetic_events(100; mission=mission)
        
        @test events.meta.headers["TELESCOP"] == uppercase(mission)
        @test events.meta.headers["INSTRUME"] == "TEST_INSTRUMENT"
        @test haskey(events.meta.headers, "NAXIS2")
        @test haskey(events.meta.headers, "TSTART")
        @test haskey(events.meta.headers, "TSTOP")
    end
end

# Test: EventList Interface Methods
let
    events = create_synthetic_events(1000)
    
    # Test length methods
    @test length(events) == 1000
    @test size(events) == (1000,)
    
    # Test accessor methods
    @test times(events) === events.times
    @test energies(events) === events.energies
    @test has_energies(events) == true
    
    # Test show methods (should not error)
    io = IOBuffer()
    show(io, MIME"text/plain"(), events)
    @test length(String(take!(io))) > 0
    
    show(io, MIME"text/plain"(), events.meta)
    @test length(String(take!(io))) > 0
end

# Test: Large Dataset Creation
let
    # Test that large datasets can be created without issues
    @time events = create_synthetic_events(50000)
    @test length(events) == 50000
    @test sizeof(times(events)) + sizeof(energies(events)) < 1e6  # Less than 1MB for 50k events
end

# Test: Memory Efficiency
let
    # Test that no unnecessary copies are made
    events1 = create_synthetic_events(1000)
    events2 = create_synthetic_events(1000)
    
    # Each should have independent data
    @test times(events1) !== times(events2)
    @test energies(events1) !== energies(events2)
    @test events1.meta.extra_columns["PI"] !== events2.meta.extra_columns["PI"]
end

# Test: Random Seed Behavior
let
    # Test that different calls produce different results (random behavior)
    events1 = create_synthetic_events(1000)
    events2 = create_synthetic_events(1000)
    
    # Should be different due to randomness
    @test times(events1) != times(events2)
    @test events1.meta.extra_columns["PI"] != events2.meta.extra_columns["PI"]
    @test events1.meta.extra_columns["DETID"] != events2.meta.extra_columns["DETID"]
    
    # But same structure
    @test length(events1) == length(events2)
    @test typeof(events1) == typeof(events2)
end

# Test: Input Validation - Error Cases
let
    # Test negative event count
    @test_throws ArgumentError create_synthetic_events(-1)
    @test_throws ArgumentError create_synthetic_events(0)
    
    # Test invalid time range
    @test_throws ArgumentError create_synthetic_events(100; time_range=(100.0, 50.0))
    @test_throws ArgumentError create_synthetic_events(100; time_range=(100.0, 100.0))
    
    # Test invalid PI range
    @test_throws ArgumentError create_synthetic_events(100; pi_range=(1000, 500))
end

# Test: FITSMetadata Structure
let
    events = create_synthetic_events(100)
    meta = events.meta
    
    # Test all required fields are present
    @test isa(meta.filepath, String)
    @test isa(meta.hdu, Int)
    @test isa(meta.energy_units, Union{String, Nothing})
    @test isa(meta.extra_columns, Dict{String, Vector})
    @test !isnothing(meta.headers)
    
    # Test metadata consistency
    @test meta.hdu == 2
    @test meta.energy_units == "ENERGY"
    @test length(meta.extra_columns) == 2  # PI and DETID
end

# Test: Summary Function
let
    events = create_synthetic_events(1000)
    summary_str = summary(events)
    
    @test isa(summary_str, String)
    @test occursin("1000 events", summary_str)
    # The time span should be close to but less than the full range (1000.0)
    actual_time_span = maximum(times(events)) - minimum(times(events))
    expected_pattern = "$(actual_time_span) time units"
    @test occursin(expected_pattern, summary_str)
    
    # Alternative: Use regex to match any floating point number
    @test occursin(r"\d+\.\d+ time units", summary_str)
    
    @test occursin("energies:", summary_str)
    @test occursin("(ENERGY)", summary_str)
    @test occursin("2 extra columns", summary_str)
end