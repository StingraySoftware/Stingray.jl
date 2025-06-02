@testset "EventList Tests" begin
    
    # Test 1: Basic EventList creation and validation
    @testset "EventList Constructor Validation" begin
        test_dir = mktempdir()
        filename = joinpath(test_dir, "test.fits")
        metadata = DictMetadata([Dict{String,Any}()])
        
        # Test valid construction with sorted data
        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        energies = [10.0, 20.0, 15.0, 25.0, 30.0]
        extra_cols = Dict{String, Vector}("DETX" => [0.1, 0.2, 0.3, 0.4, 0.5])
        
        ev = EventList{Float64}(filename, times, energies, extra_cols, metadata)
        @test ev.filename == filename
        @test ev.times == times
        @test ev.energies == energies
        @test ev.extra_columns == extra_cols
        @test ev.metadata == metadata
        
        # Test automatic sorting with unsorted data
        unsorted_times = [3.0, 1.0, 4.0, 2.0]
        unsorted_energies = [15.0, 10.0, 25.0, 20.0]
        unsorted_extra_cols = Dict{String, Vector}("DETX" => [0.3, 0.1, 0.4, 0.2])
        
        ev_unsorted = EventList{Float64}(filename, unsorted_times, unsorted_energies, unsorted_extra_cols, metadata)
        @test issorted(ev_unsorted.times)
        @test ev_unsorted.times == sort(unsorted_times)
        @test ev_unsorted.energies == [10.0, 20.0, 15.0, 25.0]  # Values should follow the sorted order
        @test ev_unsorted.extra_columns["DETX"] == [0.1, 0.2, 0.3, 0.4]  # Values should follow the sorted order
        
        # Test validation: empty times should throw
        @test_throws ArgumentError EventList{Float64}(filename, Float64[], nothing, Dict{String, Vector}(), metadata)
        
        # Test validation: mismatched energy vector length
        wrong_energies = [10.0, 20.0]  # Only 2 elements vs 4 times
        @test_throws ArgumentError EventList{Float64}(filename, unsorted_times, wrong_energies, Dict{String, Vector}(), metadata)
        
        # Test validation: mismatched extra column length
        wrong_extra = Dict{String, Vector}("DETX" => [0.1, 0.2])  # Only 2 elements vs 4 times
        @test_throws ArgumentError EventList{Float64}(filename, unsorted_times, nothing, wrong_extra, metadata)
    end
    
    # Test 2: Simplified constructors
    @testset "Simplified Constructors" begin
        test_dir = mktempdir()
        filename = joinpath(test_dir, "test.fits")
        times = [1.0, 2.0, 3.0]
        metadata = DictMetadata([Dict{String,Any}()])
        
        # Constructor with just times and metadata
        ev1 = EventList{Float64}(filename, times, metadata)
        @test ev1.filename == filename
        @test ev1.times == times
        @test isnothing(ev1.energies)
        @test isempty(ev1.extra_columns)
        @test ev1.metadata == metadata
        
        # Constructor with times, energies, and metadata
        energies = [10.0, 20.0, 30.0]
        ev2 = EventList{Float64}(filename, times, energies, metadata)
        @test ev2.filename == filename
        @test ev2.times == times
        @test ev2.energies == energies
        @test isempty(ev2.extra_columns)
        @test ev2.metadata == metadata
    end
    
    # Test 3: Accessor functions
    @testset "Accessor Functions" begin
        test_dir = mktempdir()
        filename = joinpath(test_dir, "test.fits")
        times_vec = [1.0, 2.0, 3.0]
        energies_vec = [10.0, 20.0, 30.0]
        metadata = DictMetadata([Dict{String,Any}()])
        
        ev = EventList{Float64}(filename, times_vec, energies_vec, metadata)
        
        # Test times() accessor
        @test times(ev) === ev.times
        @test times(ev) == times_vec
        
        # Test energies() accessor  
        @test energies(ev) === ev.energies
        @test energies(ev) == energies_vec
        
        # Test with nothing energies
        ev_no_energy = EventList{Float64}(filename, times_vec, metadata)
        @test isnothing(energies(ev_no_energy))
    end
    
    # Test 4: Base interface methods
    @testset "Base Interface Methods" begin
        test_dir = mktempdir()
        filename = joinpath(test_dir, "test.fits")
        times_vec = [1.0, 2.0, 3.0, 4.0]
        metadata = DictMetadata([Dict{String,Any}()])
        
        ev = EventList{Float64}(filename, times_vec, metadata)
        
        # Test length
        @test length(ev) == 4
        @test length(ev) == length(times_vec)
        
        # Test size
        @test size(ev) == (4,)
        @test size(ev) == (length(times_vec),)
        
        # Test show method
        io = IOBuffer()
        show(io, ev)
        str = String(take!(io))
        @test occursin("EventList{Float64}", str)
        @test occursin("n=4", str)
        @test occursin("no energy data", str)
        @test occursin("0 extra columns", str)
        @test occursin("file=$filename", str)
        
        # Test show with energy data
        energies_vec = [10.0, 20.0, 30.0, 40.0]
        ev_with_energy = EventList{Float64}(filename, times_vec, energies_vec, metadata)
        io2 = IOBuffer()
        show(io2, ev_with_energy)
        str2 = String(take!(io2))
        @test occursin("with energy data", str2)
    end
    
    @testset "readevents Basic Functionality" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample.fits")
        
        # Create a sample FITS file
        f = FITS(sample_file, "w")
        
        # Create primary HDU with a small array instead of empty
        write(f, [0]) # Use a single element array instead of empty
        
        # Create event table in HDU 2
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        
        # Test reading with default parameters
        data = readevents(sample_file)
        @test data.filename == sample_file
        @test data.times == times
        @test data.energies == energies
        @test eltype(data.times) == Float64
        @test eltype(data.energies) == Float64
        @test length(data.metadata.headers) >= 2
        
        # Test reading with different numeric type
        data_f32 = readevents(sample_file, T=Float32)
        @test eltype(data_f32.times) == Float32
        @test eltype(data_f32.energies) == Float32
        @test data_f32.times ≈ Float32.(times)
        @test data_f32.energies ≈ Float32.(energies)
    end
    
    @testset "readevents HDU Handling" begin
        test_dir = mktempdir()
        
        # Test with events in HDU 3 instead of default HDU 2
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
        
        # Should find events in HDU 3 via fallback mechanism
        data = readevents(sample_file)
        @test data.times == times
        @test data.energies == energies
        
        # Test specifying specific HDU
        data_hdu3 = readevents(sample_file, event_hdu=3)
        @test data_hdu3.times == times
        @test data_hdu3.energies == energies
    end
    
    @testset "readevents Alternative Energy Columns" begin
        test_dir = mktempdir()
        
        # Test with PI column
        pi_file = joinpath(test_dir, "pi_sample.fits")
        f = FITS(pi_file, "w")
        write(f, [0])  # Non-empty primary HDU
        
        times = Float64[1.0, 2.0, 3.0]
        pi_values = Float64[100.0, 200.0, 300.0]
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["PI"] = pi_values
        write(f, table)
        close(f)
        
        data = readevents(pi_file)
        @test data.times == times
        @test data.energies == pi_values
        
        # Test with PHA column
        pha_file = joinpath(test_dir, "pha_sample.fits")
        f = FITS(pha_file, "w")
        write(f, [0])  # Non-empty primary HDU
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["PHA"] = pi_values
        write(f, table)
        close(f)
        
        data_pha = readevents(pha_file)
        @test data_pha.times == times
        @test data_pha.energies == pi_values
    end
    
    @testset "readevents Missing Columns" begin
        test_dir = mktempdir()
        
        # File with only TIME column
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
        
        # File with no TIME column should throw error
        no_time_file = joinpath(test_dir, "no_time.fits")
        f = FITS(no_time_file, "w")
        write(f, [0])  # Non-empty primary HDU
        
        table = Dict{String,Array}()
        table["ENERGY"] = Float64[10.0, 20.0, 30.0]
        write(f, table)
        close(f)
        
        @test_throws ArgumentError readevents(no_time_file)
    end
    
    @testset "readevents Extra Columns" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "extra_cols.fits")
        
        f = FITS(sample_file, "w")
        write(f, [0])  # Non-empty primary HDU
        
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        sectors = Int64[1, 2, 1]
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        table["SECTOR"] = sectors
        write(f, table)
        close(f)
        
        # Test reading with sector column specified
        data = readevents(sample_file, sector_column="SECTOR")
        @test data.times == times
        @test data.energies == energies
        @test haskey(data.extra_columns, "SECTOR")
        @test data.extra_columns["SECTOR"] == sectors
    end
    
    # Test 10: Error handling
   @testset "Error Handling" begin
        test_dir = mktempdir()
        
        # Test non-existent file
        @test_throws CFITSIO.CFITSIOError readevents("non_existent_file.fits")
        
        # Test invalid FITS file
        invalid_file = joinpath(test_dir, "invalid.fits")
        open(invalid_file, "w") do io
            write(io, "This is not a FITS file")
        end
        @test_throws Exception readevents(invalid_file)
        
        # Test with non-table HDU specified
        sample_file = joinpath(test_dir, "image_hdu.fits")
        f = FITS(sample_file, "w")
        
        # Create a valid primary HDU with a small image
        primary_data = reshape([1.0], 1, 1)  # 1x1 image instead of empty array
        write(f, primary_data)
        
        # Create an image HDU
        image_data = reshape(collect(1:100), 10, 10)
        write(f, image_data)
        close(f)
        
        @test_throws ArgumentError readevents(sample_file, event_hdu=2)
    end
    
    @testset "Case Insensitive Column Names" begin
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
    
    @testset "Integration Test" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "realistic.fits")
        
        # Create more realistic test data
        f = FITS(sample_file, "w")
        
        # Primary HDU with proper header
        primary_data = reshape([1.0], 1, 1)  # Use 1x1 image
        header_keys = ["TELESCOP", "INSTRUME"]
        header_values = ["TEST_SAT", "TEST_DET"]
        header_comments = ["Test telescope", "Test detector"]
        primary_hdr = FITSHeader(header_keys, header_values, header_comments)
        write(f, primary_data; header=primary_hdr)
        
        # Event data with realistic values
        n_events = 1000
        times = sort(rand(n_events) * 1000.0)  # 1000 seconds of data
        energies = rand(n_events) * 10.0 .+ 0.5  # 0.5-10.5 keV
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        
        # Create event HDU header
        event_header_keys = ["EXTNAME", "TELESCOP"]
        event_header_values = ["EVENTS", "TEST_SAT"]
        event_header_comments = ["Extension name", "Test telescope"]
        event_hdr = FITSHeader(event_header_keys, event_header_values, event_header_comments)
        write(f, table; header=event_hdr)
        close(f)
        
        # Test reading
        data = readevents(sample_file)
        
        @test length(data.times) == n_events
        @test length(data.energies) == n_events
        @test issorted(data.times)
        @test minimum(data.energies) >= 0.5
        @test maximum(data.energies) <= 10.5
        @test length(data.metadata.headers) == 2
        @test data.metadata.headers[1]["TELESCOP"] == "TEST_SAT"
        
        # Check if EXTNAME exists in the second header
        if haskey(data.metadata.headers[2], "EXTNAME")
            @test data.metadata.headers[2]["EXTNAME"] == "EVENTS"
        else
            @test data.metadata.headers[2]["TELESCOP"] == "TEST_SAT"
        end
    end
end