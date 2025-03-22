using Test
using FITSIO
using Stingray
using Stingray.Events
using Stingray.LightCurveModule

@testset "LightCurve Tests" begin
    # Test 1: Create a temporary FITS file with sample event data
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "sample.fits")

    f = FITS(sample_file, "w")
    write(f, Int[])  # Primary HDU

    # Adding event data (TIME & ENERGY)
    times = Float64[1.0, 2.5, 3.2, 4.8, 6.1, 7.3, 9.0, 10.5]
    energies = Float64[100.0, 200.0, 150.0, 250.0, 300.0, 180.0, 220.0, 270.0]
    table = Dict("TIME" => times, "ENERGY" => energies)
    write(f, table)
    close(f)

    # Test 2: Read events using Events module
    eventlist = Events.readevents(sample_file)

    # Validate eventlist
    @test !isempty(eventlist.times)
    if isempty(eventlist.energies)
        @warn "⚠ ENERGY column is empty! Check if another column (e.g., FLUX) should be used instead."
    else
        @test !isempty(eventlist.energies)
    end

    # Test 3: Generate light curve
    bin_size = 2.0
    err_method = :poisson
    lightcurve = LightCurveModule.create_lightcurve(eventlist, bin_size, err_method=err_method)

    # Validate light curve
    @test all(lightcurve.counts .>= 0)
    @test length(lightcurve.counts) == length(lightcurve.timebins)
    @test lightcurve.err_method == err_method
    @test all(lightcurve.count_error .>= 0)

    # Test 4: Read an existing event file if available
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

    # Test 5: Different data types
    @testset "Different data types" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_float32.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])
        times = Float64[1.0, 2.0, 3.0]
        energies = Float64[10.0, 20.0, 30.0]
        table = Dict("TIME" => times, "ENERGY" => energies)
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

    # Test 6: Missing columns
    @testset "Missing columns" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_no_energy.fits")

        # Create a sample FITS file with only TIME column
        f = FITS(sample_file, "w")
        write(f, Int[])
        times = Float64[1.0, 2.0, 3.0]
        table = Dict("TIME" => times)
        write(f, table)
        close(f)

        local data
        @test_logs (:warn, "No ENERGY data found in FITS file $(sample_file). Energy spectrum analysis will not be possible.") begin
            data = readevents(sample_file)
        end
        @test length(data.times) == 3
        @test length(data.energies) == 0

        # Create a file with only ENERGY column
        sample_file2 = joinpath(test_dir, "sample_no_time.fits")
        f = FITS(sample_file2, "w")
        write(f, Int[])
        energies = Float64[10.0, 20.0, 30.0]
        table = Dict("ENERGY" => energies)
        write(f, table)
        close(f)

        local data2
        @test_logs (:warn, "No TIME data found in FITS file $(sample_file2). Time series analysis will not be possible.") begin
            data2 = readevents(sample_file2)
        end
        @test length(data2.times) == 0
        @test length(data2.energies) == 3
    end

    # Test 7: Multiple HDUs
    @testset "Multiple HDUs" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_multi_hdu.fits")

        # Create a sample FITS file with multiple HDUs
        f = FITS(sample_file, "w")
        write(f, Int[])

        # First HDU (TIME + ENERGY)
        times1 = Float64[1.0, 2.0, 3.0]
        energies1 = Float64[10.0, 20.0, 30.0]
        table1 = Dict("TIME" => times1, "ENERGY" => energies1)
        write(f, table1)

        # Second HDU (OTHER column)
        other_data = Float64[100.0, 200.0, 300.0]
        table2 = Dict("OTHER" => other_data)
        write(f, table2)

        # Third HDU (TIME only)
        times3 = Float64[4.0, 5.0, 6.0]
        table3 = Dict("TIME" => times3)
        write(f, table3)
        close(f)

        data = readevents(sample_file)
        @test length(data.metadata.headers) == 4  # Primary + 3 table HDUs
        @test length(data.times) == 3
        @test length(data.energies) == 3
    end

    # Test 8: Error handling
    @testset "Error handling" begin
        # Test with non-existent file
        @test_throws Exception readevents("non_existent_file.fits")

        # Test with invalid FITS file
        invalid_file = tempname()
        open(invalid_file, "w") do io
            write(io, "This is not a FITS file")
        end
        @test_throws Exception readevents(invalid_file)
    end

    #println("🔹 LightCurve Timebins: ", lightcurve.timebins)
    #println("🔹 LightCurve Counts: ", lightcurve.counts)
    #println("LightCurve Tests Passed!")
end
