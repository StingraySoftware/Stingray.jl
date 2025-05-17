@testset "Complete LightCurve Tests" begin
    @testset "Basic Light Curve Creation" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])
        
        # Create test data
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        
        data = readevents(sample_file)
        
        # Create light curve
        lc = create_lightcurve(data, 1.0)
        
        # Calculate expected bins
        expected_bins = Int(ceil((maximum(times) - minimum(times))/1.0))
        
        # Test structure
        @test length(lc.timebins) == expected_bins
        @test length(lc.counts) == expected_bins
        @test length(lc.bin_edges) == expected_bins + 1
        
        # Test bin centers
        @test lc.timebins[1] ≈ 1.5
        @test lc.timebins[end] ≈ 4.5
        
        # Test counts
        expected_counts = fill(1, expected_bins)
        @test all(lc.counts .== expected_counts)
        
        # Test errors
        @test all(lc.count_error .≈ sqrt.(Float64.(expected_counts)))
        
        # Test metadata and properties
        @test lc.err_method === :poisson
        @test length(lc) == expected_bins
        @test size(lc) == (expected_bins,)
        @test lc[1] == (1.5, 1)
    end

    @testset "Time Range and Binning" begin
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
        events = EventList{Float64}("test.fits", times, energies, 
            DictMetadata([Dict{String,Any}()]))

        # Test specific time range
        lc = create_lightcurve(events, 1.0, tstart=2.0, tstop=4.0)
        expected_bins = Int(ceil((4.0 - 2.0)/1.0))
        @test length(lc.timebins) == expected_bins
        @test lc.metadata.time_range == (2.0, 4.0)
        @test all(2.0 .<= lc.bin_edges .<= 4.0)
        @test sum(lc.counts) == 2

        # Test equal start and stop times
        lc_equal = create_lightcurve(events, 1.0, tstart=2.0, tstop=2.0)
        @test length(lc_equal.counts) == 1
        @test lc_equal.metadata.time_range[2] == lc_equal.metadata.time_range[1] + 1.0

        # Test bin edges
        lc_edges = create_lightcurve(events, 2.0)
        @test lc_edges.bin_edges[end] >= maximum(times)
    end

    @testset "Filtering" begin
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[1.0, 2.0, 5.0, 8.0, 10.0]
        events = EventList{Float64}("test.fits", times, energies,
            DictMetadata([Dict{String,Any}()]))

        # Test energy filtering
        energy_filter = Dict{Symbol,Any}(:energy => (4.0, 9.0))
        lc = create_lightcurve(events, 1.0, filters=energy_filter)
        @test sum(lc.counts) == 2
        @test haskey(lc.metadata.extra, "filtered_nevents")
        @test lc.metadata.extra["filtered_nevents"] == 2

        # Test empty filter result
        empty_filter = Dict{Symbol,Any}(:energy => (100.0, 200.0))
        lc_empty = create_lightcurve(events, 1.0, filters=empty_filter)
        @test all(lc_empty.counts .== 0)
        @test haskey(lc_empty.metadata.extra, "warning")
        @test lc_empty.metadata.extra["warning"] == "No events remain after filtering"
    end

    @testset "Error Methods" begin
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        events = EventList{Float64}("test.fits", times, energies,
            DictMetadata([Dict{String,Any}()]))

        # Test Poisson errors
        lc_poisson = create_lightcurve(events, 1.0)
        @test all(lc_poisson.count_error .≈ sqrt.(lc_poisson.counts))

        # Test Gaussian errors
        lc_gaussian = create_lightcurve(events, 1.0, err_method=:gaussian)
        @test all(lc_gaussian.count_error .≈ sqrt.(lc_gaussian.counts .+ 1))

        # Test invalid error method
        @test_throws ArgumentError create_lightcurve(events, 1.0, 
            err_method=:invalid)
    end

    @testset "Properties and Metadata" begin
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        headers = [Dict{String,Any}(
            "TELESCOP" => "TEST",
            "INSTRUME" => "INST",
            "OBJECT" => "SRC",
            "MJDREF" => 58000.0
        )]
        events = EventList{Float64}("test.fits", times, energies,
            DictMetadata(headers))

        lc = create_lightcurve(events, 1.0)
        
        # Test metadata
        @test lc.metadata.telescope == "TEST"
        @test lc.metadata.instrument == "INST"
        @test lc.metadata.object == "SRC"
        @test lc.metadata.mjdref == 58000.0
        @test haskey(lc.metadata.extra, "filtered_nevents")
        @test haskey(lc.metadata.extra, "total_nevents")

        # Test properties
        @test !isempty(lc.properties)
        energy_prop = first(filter(p -> p.name == :mean_energy, lc.properties))
        @test energy_prop.unit == "keV"
        @test length(energy_prop.values) == length(lc.counts)
    end

    @testset "Rebinning" begin
        # Create test data with evenly spaced events
        times = Float64[1.0, 1.5, 2.0, 2.5, 3.0]
        energies = Float64[10.0, 15.0, 20.0, 25.0, 30.0]
        events = EventList{Float64}("test.fits", times, energies,
            DictMetadata([Dict{String,Any}()]))
    
        # Create initial light curve with 0.5 bin size
        lc = create_lightcurve(events, 0.5)
        
        # Calculate expected number of bins after rebinning
        time_range = lc.metadata.time_range[2] - lc.metadata.time_range[1]
        expected_bins = Int(ceil(time_range))  # For 1.0 binsize
        
        # Test rebinning
        lc_rebinned = rebin(lc, 1.0)
        @test length(lc_rebinned.counts) == expected_bins
        @test sum(lc_rebinned.counts) == sum(lc.counts)
        @test all(lc_rebinned.exposure .== 1.0)
    
        # Test property preservation
        @test length(lc_rebinned.properties) == length(lc.properties)
        if !isempty(lc.properties)
            orig_prop = first(lc.properties)
            rebin_prop = first(lc_rebinned.properties)
            @test orig_prop.name == rebin_prop.name
            @test orig_prop.unit == rebin_prop.unit
        end
    
        # Test metadata
        @test haskey(lc_rebinned.metadata.extra, "original_binsize")
        @test lc_rebinned.metadata.extra["original_binsize"] == 0.5
    
        # Test invalid rebin size
        @test_throws ArgumentError rebin(lc, 0.1)
    end
    
    @testset "Edge Cases" begin
        # Test empty event list
        empty_events = EventList{Float64}("test.fits", Float64[], Float64[],
            DictMetadata([Dict{String,Any}()]))
        @test_throws ArgumentError create_lightcurve(empty_events, 1.0)
    
        # Test single event
        # Place event exactly at bin center to ensure it's counted
        times = Float64[2.5]  # Place at 2.5 to ensure it falls in a bin center
        energies = Float64[10.0]
        single_event = EventList{Float64}("test.fits", times, energies,
            DictMetadata([Dict{String,Any}()]))
        
        lc_single = create_lightcurve(single_event, 1.0)
        
        # Calculate expected bin for the event
        start_time = floor(minimum(times))
        bin_idx = Int(floor((times[1] - start_time) / 1.0)) + 1
        expected_counts = zeros(Int, length(lc_single.counts))
        if 1 <= bin_idx <= length(expected_counts)
            expected_counts[bin_idx] = 1
        end
        
        @test lc_single.counts == expected_counts
        @test sum(lc_single.counts) == 1
    
        # Test invalid bin sizes
        events = EventList{Float64}("test.fits", [1.0, 2.0], [10.0, 20.0],
            DictMetadata([Dict{String,Any}()]))
        @test_throws ArgumentError create_lightcurve(events, 0.0)
        @test_throws ArgumentError create_lightcurve(events, -1.0)
    
        # Test complete filtering
        lc_filtered = create_lightcurve(events, 1.0,
            filters=Dict{Symbol,Any}(:energy => (100.0, 200.0)))
        @test all(lc_filtered.counts .== 0)
        @test haskey(lc_filtered.metadata.extra, "warning")
    end

    @testset "Type Stability" begin
        for T in [Float32, Float64]
            times = T[1.0, 2.0, 3.0]
            energies = T[10.0, 20.0, 30.0]
            events = EventList{T}("test.fits", times, energies,
                DictMetadata([Dict{String,Any}()]))

            # Test creation
            lc = create_lightcurve(events, T(1.0))
            @test eltype(lc.timebins) === T
            @test eltype(lc.bin_edges) === T
            @test eltype(lc.count_error) === T
            @test eltype(lc.exposure) === T

            # Test rebinning
            lc_rebinned = rebin(lc, T(2.0))
            @test eltype(lc_rebinned.timebins) === T
            @test eltype(lc_rebinned.bin_edges) === T
            @test eltype(lc_rebinned.count_error) === T
            @test eltype(lc_rebinned.exposure) === T
        end
    end

    @testset "Array Interface" begin
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        events = EventList{Float64}("test.fits", times, energies,
            DictMetadata([Dict{String,Any}()]))
        
        lc = create_lightcurve(events, 1.0)
        
        @test length(lc) == length(lc.counts)
        @test size(lc) == (length(lc.counts),)
        @test lc[1] == (lc.timebins[1], lc.counts[1])
    end
end