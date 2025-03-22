@testset "EventList Tests" begin
    # Test 1: Create a sample FITS file for testing
    @testset "Sample FITS file creation" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])  # Empty primary array
        # Create a binary table HDU with TIME and ENERGY columns
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
        # Add a binary table extension
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)

        @test isfile(sample_file)

        # Test reading the sample file
        data = readevents(sample_file)
        @test data.filename == sample_file
        @test length(data.times) == 5
        @test length(data.energies) == 5
        @test eltype(data.times) == Float64
        @test eltype(data.energies) == Float64
    end

    # Test 2: Test with different data types
    @testset "Different data types" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_float32.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])
        # Create data
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        # Test with Float32
        data_f32 = readevents(sample_file, T = Float32)
        @test eltype(data_f32.times) == Float32
        @test eltype(data_f32.energies) == Float32
        # Test with Int64
        data_i64 = readevents(sample_file, T = Int64)
        @test eltype(data_i64.times) == Int64
        @test eltype(data_i64.energies) == Int64
    end

    # Test 3: Test with missing columns
    @testset "Missing columns" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_no_energy.fits")
        # Create a sample FITS file with only TIME column
        f = FITS(sample_file, "w")
        write(f, Int[])
        times = Float64[1.0, 2.0, 3.0]
        table = Dict{String,Array}()
        table["TIME"] = times
        write(f, table)
        close(f)
        local data
        @test_logs (:warn, "No ENERGY data found in FITS file $(sample_file). Energy spectrum analysis will not be possible.") begin
            data = readevents(sample_file)
        end
        @test length(data.times) == 3
        @test length(data.energies) == 0
        #create a file with only ENERGY column
        sample_file2 = joinpath(test_dir, "sample_no_time.fits")
        f = FITS(sample_file2, "w")
        write(f, Int[])  # Empty primary array
        energies = Float64[10.0, 20.0, 30.0]
        table = Dict{String,Array}()
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        local data2
        @test_logs (:warn, "No TIME data found in FITS file $(sample_file2). Time series analysis will not be possible.") begin
            data2 = readevents(sample_file2)
        end
        @test length(data2.times) == 0  # No TIME column
        @test length(data2.energies) == 3
    end

    # Test 4: Test with multiple HDUs
    @testset "Multiple HDUs" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_multi_hdu.fits")
        # Create a sample FITS file with multiple HDUs
        f = FITS(sample_file, "w")
        write(f, Int[])
        times1 = Float64[1.0, 2.0, 3.0]
        energies1 = Float64[10.0, 20.0, 30.0]
        table1 = Dict{String,Array}()
        table1["TIME"] = times1
        table1["ENERGY"] = energies1
        write(f, table1)
        # Second table HDU (with OTHER column)
        other_data = Float64[100.0, 200.0, 300.0]
        table2 = Dict{String,Array}()
        table2["OTHER"] = other_data
        write(f, table2)
        # Third table HDU (with TIME only)
        times3 = Float64[4.0, 5.0, 6.0]
        table3 = Dict{String,Array}()
        table3["TIME"] = times3
        write(f, table3)
        close(f)
        data = readevents(sample_file)
        @test length(data.metadata.headers) == 4  # Primary + 3 table HDUs
        # Should read the first HDU with both TIME and ENERGY
        @test length(data.times) == 3
        @test length(data.energies) == 3
    end

    # Test 5: Test with monol_testA.evt
    @testset "test monol_testA.evt" begin
        test_filepath = joinpath("data", "monol_testA.evt")
        if isfile(test_filepath)
            @testset "monol_testA.evt" begin
                old_logger = global_logger(ConsoleLogger(stderr, Logging.Error))
                try
                    data = readevents(test_filepath)
                    @test data.filename == test_filepath
                    @test length(data.metadata.headers) > 0
                finally
                    global_logger(old_logger)
                end
            end
        else
            @info "Test file '$(test_filepath)' not found. Skipping this test."
        end
    end

    @testset "Error handling" begin
        # Test with non-existent file - using a more generic approach
        @test_throws Exception readevents("non_existent_file.fits")

        # Test with invalid FITS file
        invalid_file = tempname()
        open(invalid_file, "w") do io
            write(io, "This is not a FITS file")
        end
        @test_throws Exception readevents(invalid_file)
    end
end
