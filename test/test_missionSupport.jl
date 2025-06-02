"""
    create_synthetic_events(n_events::Int=1000; mission::String="nustar", 
                           time_range::Tuple{Float64,Float64}=(0.0, 1000.0),
                           pi_range::Tuple{Int,Int}=(1, 4096),
                           T::Type=Float64) -> EventList{T}

Create synthetic X-ray event data for testing purposes.

This function generates realistic synthetic X-ray event data with proper time ordering,
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
- `EventList{T}`: Synthetic event list with times, energies, and metadata

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
    filename = "synthetic_$(mission).fits"
    
    # Create metadata headers
    main_header = Dict{String, Any}(
        "TELESCOP" => uppercase(mission),
        "INSTRUME" => "TEST_INSTRUMENT",
        "OBSERVER" => "SYNTHETIC_DATA",
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
        main_header["DETNAM"] = "TEST_DET"
        main_header["PIFLTCOR"] = "T"
    elseif lowercase(mission) == "xmm"
        main_header["FILTER"] = "Medium"
        main_header["SUBMODE"] = "PrimeFullWindow"
    elseif lowercase(mission) == "nicer"
        main_header["DETNAM"] = "TEST_NICER"
        main_header["FILTFILE"] = "NONE"
    elseif lowercase(mission) in ["chandra", "axaf"]
        main_header["DETNAM"] = "ACIS-S"
        main_header["GRATING"] = "NONE"
    elseif lowercase(mission) == "xte"
        main_header["DETNAM"] = "PCA"
        main_header["LAYERS"] = "ALL"
    elseif lowercase(mission) == "ixpe"
        main_header["DETNAM"] = "GPD"
        main_header["POLMODE"] = "ON"
    end
    
    # Create additional headers (empty for now, but maintaining structure)
    headers = [main_header]
    
    # Create DictMetadata object
    metadata = DictMetadata(headers)
    
    # Create and return EventList
    return EventList{T}(filename, times, energies, extra_columns, metadata)
end

@testset verbose=true "Synthetic Events Tests" begin
    
    @testset "Basic Synthetic Event Creation" begin
        @testset "Default Parameters" begin
            events = create_synthetic_events()
            
            # Test basic structure
            @test isa(events, EventList)
            @test length(events.times) == 1000  # Default n_events
            @test length(events.energies) == 1000
            @test events.filename == "synthetic_nustar.fits"
            
            # Test times are sorted and in range
            @test issorted(events.times)
            @test all(0.0 .<= events.times .<= 1000.0)
            
            # Test energies are positive (after calibration)
            @test all(events.energies .> 0)
            
            # Test extra columns
            @test haskey(events.extra_columns, "PI")
            @test haskey(events.extra_columns, "DETID")
            @test length(events.extra_columns["PI"]) == 1000
            @test length(events.extra_columns["DETID"]) == 1000
            
            # Test PI channels are in expected range
            @test all(1 .<= events.extra_columns["PI"] .<= 4096)
            
            # Test detector IDs are in expected range
            @test all(0 .<= events.extra_columns["DETID"] .<= 3)
            
            # Test metadata structure
            @test isa(events.metadata, DictMetadata)
            @test length(events.metadata.headers) >= 1
            
            # Test metadata content (assuming first header contains main info)
            main_header = events.metadata.headers[1]
            @test main_header["TELESCOP"] == "NUSTAR"
            @test main_header["INSTRUME"] == "TEST_INSTRUMENT"
            @test main_header["NAXIS2"] == 1000
            @test main_header["TSTART"] == 0.0
            @test main_header["TSTOP"] == 1000.0
        end
        
        @testset "Custom Parameters" begin
            n_events = 500
            mission = "xmm"
            time_range = (100.0, 200.0)
            pi_range = (50, 1000)
            
            events = create_synthetic_events(n_events; 
                                           mission=mission, 
                                           time_range=time_range, 
                                           pi_range=pi_range)
            
            # Test custom parameters are respected
            @test length(events.times) == n_events
            @test length(events.energies) == n_events
            @test events.filename == "synthetic_xmm.fits"
            
            # Test time range
            @test all(time_range[1] .<= events.times .<= time_range[2])
            main_header = events.metadata.headers[1]
            @test main_header["TSTART"] == time_range[1]
            @test main_header["TSTOP"] == time_range[2]
            
            # Test PI range
            @test all(pi_range[1] .<= events.extra_columns["PI"] .<= pi_range[2])
            
            # Test mission-specific metadata
            @test main_header["TELESCOP"] == "XMM"
        end
        
        @testset "Different Missions" begin
            missions = ["nustar", "xmm", "nicer", "ixpe", "axaf", "chandra", "xte"]
            
            for mission in missions
                events = create_synthetic_events(100; mission=mission)
                
                @test events.filename == "synthetic_$(mission).fits"
                main_header = events.metadata.headers[1]
                @test main_header["TELESCOP"] == uppercase(mission)
                @test length(events.times) == 100
                @test length(events.energies) == 100
                
                # Test that calibration was applied correctly
                ms = get_mission_support(mission)
                expected_energies = apply_calibration(ms, events.extra_columns["PI"])
                @test events.energies ≈ expected_energies
            end
        end
    end
    
    @testset "Data Quality and Consistency" begin
        @testset "Time Ordering" begin
            for _ in 1:10  # Test multiple times due to randomness
                events = create_synthetic_events(100)
                @test issorted(events.times)
                
                # Test no duplicate times (very unlikely but possible)
                @test length(unique(events.times)) >= 95  # Allow some duplicates due to floating point
            end
        end
        
        @testset "Energy Calibration Consistency" begin
            missions = ["nustar", "xmm", "nicer", "ixpe"]
            
            for mission in missions
                events = create_synthetic_events(200; mission=mission)
                
                # Manually verify calibration
                ms = get_mission_support(mission)
                expected_energies = apply_calibration(ms, events.extra_columns["PI"])
                @test events.energies ≈ expected_energies
                
                # Test energy ranges are reasonable for each mission
                if mission == "nustar"
                    # NuSTAR: pi * 0.04 + 1.62, PI range 1-4096
                    @test minimum(events.energies) >= 1.66  # 1*0.04 + 1.62
                    @test maximum(events.energies) <= 165.46  # 4096*0.04 + 1.62
                elseif mission == "xmm"
                    # XMM: pi * 0.001, PI range 1-4096
                    @test minimum(events.energies) >= 0.001
                    @test maximum(events.energies) <= 4.096
                elseif mission == "nicer"
                    # NICER: pi * 0.01, PI range 1-4096
                    @test minimum(events.energies) >= 0.01
                    @test maximum(events.energies) <= 40.96
                elseif mission == "ixpe"
                    # IXPE: pi / 375 * 15, PI range 1-4096
                    @test minimum(events.energies) >= 15.0/375  # ≈ 0.04
                    @test maximum(events.energies) <= 4096*15.0/375  # ≈ 163.84
                end
            end
        end
        
        @testset "Statistical Properties" begin
            events = create_synthetic_events(10000)  # Large sample for statistics
            
            # Test time distribution (should be roughly uniform)
            time_hist = fit(Histogram, events.times, 0:100:1000)
            counts = time_hist.weights
            # Expect roughly equal counts in each bin (within statistical fluctuation)
            expected_count = 10000 / 10  # 10 bins
            @test all(abs.(counts .- expected_count) .< 3 * sqrt(expected_count))  # 3-sigma test
            
            # Test PI distribution (should be roughly uniform over discrete range)
            pi_hist = fit(Histogram, events.extra_columns["PI"], 1:100:4096)
            pi_counts = pi_hist.weights
            # More lenient test due to discrete uniform distribution
            @test std(pi_counts) / mean(pi_counts) < 0.2  # Coefficient of variation < 20%
            
            # Test detector ID distribution
            detid_counts = [count(==(i), events.extra_columns["DETID"]) for i in 0:3]
            @test all(abs.(detid_counts .- 2500) .< 3 * sqrt(2500))  # Each detector ~2500 events
        end
    end
    
    @testset "Edge Cases and Error Handling" begin
        @testset "Small Event Counts" begin
            for n in [1, 2, 5, 10]
                events = create_synthetic_events(n)
                @test length(events.times) == n
                @test length(events.energies) == n
                @test length(events.extra_columns["PI"]) == n
                @test length(events.extra_columns["DETID"]) == n
                @test issorted(events.times)
            end
        end
        
        @testset "Large Event Counts" begin
            events = create_synthetic_events(100000)
            @test length(events.times) == 100000
            @test length(events.energies) == 100000
            @test issorted(events.times)
            main_header = events.metadata.headers[1]
            @test main_header["NAXIS2"] == 100000
        end
        
        @testset "Extreme Time Ranges" begin
            # Very short time range
            events = create_synthetic_events(100; time_range=(0.0, 0.1))
            @test all(0.0 .<= events.times .<= 0.1)
            main_header = events.metadata.headers[1]
            @test main_header["TSTART"] == 0.0
            @test main_header["TSTOP"] == 0.1
            
            # Very long time range
            events = create_synthetic_events(100; time_range=(0.0, 1e6))
            @test all(0.0 .<= events.times .<= 1e6)
            main_header = events.metadata.headers[1]
            @test main_header["TSTART"] == 0.0
            @test main_header["TSTOP"] == 1e6
            
            # Negative time range
            events = create_synthetic_events(100; time_range=(-1000.0, -500.0))
            @test all(-1000.0 .<= events.times .<= -500.0)
            @test issorted(events.times)
        end
        
        @testset "Extreme PI Ranges" begin
            # Small PI range
            events = create_synthetic_events(100; pi_range=(100, 110))
            @test all(100 .<= events.extra_columns["PI"] .<= 110)
            
            # Single PI value
            events = create_synthetic_events(100; pi_range=(500, 500))
            @test all(events.extra_columns["PI"] .== 500)
            @test all(events.energies .== events.energies[1])  # All same energy
            
            # Large PI range
            events = create_synthetic_events(100; pi_range=(1, 10000))
            @test all(1 .<= events.extra_columns["PI"] .<= 10000)
        end
        
        @testset "Unknown Mission" begin
            # Should still work but with warning
            @test_logs (:warn, r"Mission unknown_mission not recognized") begin
                events = create_synthetic_events(100; mission="unknown_mission")
                @test events.filename == "synthetic_unknown_mission.fits"
                main_header = events.metadata.headers[1]
                @test main_header["TELESCOP"] == "UNKNOWN_MISSION"
                @test length(events.times) == 100
                # Should use identity calibration
                @test events.energies == Float64.(events.extra_columns["PI"])
            end
        end
    end
    
    @testset "Data Integrity" begin
        @testset "No Missing Data" begin
            events = create_synthetic_events(1000)
            
            # Check no NaN or missing values
            @test all(isfinite.(events.times))
            @test all(isfinite.(events.energies))
            @test all(isfinite.(events.extra_columns["PI"]))
            @test all(isfinite.(events.extra_columns["DETID"]))
            
            # Check no negative energies (after calibration)
            @test all(events.energies .>= 0)
        end
        
        @testset "Correct Data Types" begin
            events = create_synthetic_events(100)
            
            @test eltype(events.times) == Float64
            @test eltype(events.energies) == Float64
            @test eltype(events.extra_columns["PI"]) <: Integer
            @test eltype(events.extra_columns["DETID"]) <: Integer
            @test isa(events.metadata, DictMetadata)
            @test isa(events.filename, String)
        end
        
        @testset "Array Length Consistency" begin
            for n in [10, 100, 1000, 5000]
                events = create_synthetic_events(n)
                
                @test length(events.times) == n
                @test length(events.energies) == n
                @test length(events.extra_columns["PI"]) == n
                @test length(events.extra_columns["DETID"]) == n
                main_header = events.metadata.headers[1]
                @test main_header["NAXIS2"] == n
            end
        end
    end
    
    @testset "Mission-Specific Behavior" begin
        @testset "Mission Name Handling" begin
            # Test case sensitivity
            missions = ["NUSTAR", "nustar", "NuSTAR", "NuStar"]
            for mission in missions
                events = create_synthetic_events(100; mission=mission)
                main_header = events.metadata.headers[1]
                @test main_header["TELESCOP"] == "NUSTAR"
                @test events.filename == "synthetic_$(mission).fits"
            end
        end
        
        @testset "Calibration Differences" begin
            # Same PI values should give different energies for different missions
            test_pi_range = (1000, 1000)  # Fixed PI value
            
            results = Dict{String, Float64}()
            for mission in ["nustar", "xmm", "nicer", "ixpe"]
                events = create_synthetic_events(100; mission=mission, pi_range=test_pi_range)
                results[mission] = events.energies[1]  # All should be same since PI is fixed
            end
            
            # Different missions should give different energies
            missions = collect(keys(results))
            for i in 1:length(missions)
                for j in (i+1):length(missions)
                    @test results[missions[i]] ≠ results[missions[j]]
                end
            end
        end
        
        @testset "Metadata Consistency" begin
            missions = ["nustar", "xmm", "nicer", "ixpe", "axaf", "chandra"]
            
            for mission in missions
                events = create_synthetic_events(100; mission=mission)
                
                @test events.metadata.headers[1]["TELESCOP"] == uppercase(mission)
                @test events.metadata.headers[1]["INSTRUME"] == "TEST_INSTRUMENT"
                @test haskey(events.metadata.headers[1], "NAXIS2")
                @test haskey(events.metadata.headers[1], "TSTART")
                @test haskey(events.metadata.headers[1], "TSTOP")
            end
        end
    end
    
    @testset "Performance and Memory" begin
        @testset "Large Dataset Creation" begin
            # Test that large datasets can be created without issues
            @time events = create_synthetic_events(50000)
            @test length(events.times) == 50000
            @test sizeof(events.times) + sizeof(events.energies) < 1e6  # Less than 1MB for 50k events
        end
        
        @testset "Memory Efficiency" begin
            # Test that no unnecessary copies are made
            events1 = create_synthetic_events(1000)
            events2 = create_synthetic_events(1000)
            
            # Each should have independent data
            @test events1.times !== events2.times
            @test events1.energies !== events2.energies
            @test events1.extra_columns["PI"] !== events2.extra_columns["PI"]
        end
    end
    
    @testset "Reproducibility" begin
        @testset "Random Seed Behavior" begin
            # Test that different calls produce different results (random behavior)
            events1 = create_synthetic_events(1000)
            events2 = create_synthetic_events(1000)
            
            # Should be different due to randomness
            @test events1.times != events2.times
            @test events1.extra_columns["PI"] != events2.extra_columns["PI"]
            @test events1.extra_columns["DETID"] != events2.extra_columns["DETID"]
            
            # But same structure
            @test length(events1.times) == length(events2.times)
            @test typeof(events1) == typeof(events2)
        end
    end
end