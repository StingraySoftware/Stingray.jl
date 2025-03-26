@testset "EventList Tests" begin
    @testset "Basic functionality" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "basic.fits")
        
        FITS(sample_file, "w") do f
            write(f, Float32[])  # Empty primary HDU
            
            # Table HDU with TIME and ENERGY
            data = Dict{String,Vector{Float64}}(
                "TIME" => Float64[1:10...],
                "ENERGY" => Float64[11:20...]
            )
            write(f, data)
            
            # Additional HDU that should be ignored
            other_data = Dict{String,Vector{Float64}}(
                "RATE" => Float64[21:30...]
            )
            write(f, data)  # Write the same data again
        end

        event_list = readevents(sample_file)
        @test event_list.filename == sample_file
        @test length(event_list.times) == 10
        @test length(event_list.energies) == 10
        @test event_list.times == collect(1:10)
        @test event_list.energies == collect(11:20)
        @test length(event_list.metadata.headers) == 2  # Primary + event HDU
        
        rm(test_dir, recursive=true, force=true)
    end

    @testset "Different data types" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "datatypes.fits")
        
        FITS(sample_file, "w") do f
            write(f, Float32[])
            data = Dict{String,Vector{Float64}}(
                "TIME" => Float64[1.0, 2.0, 3.0],
                "ENERGY" => Float64[10.0, 20.0, 30.0]
            )
            write(f, data)
        end

        # Test Float32 conversion
        data_f32 = readevents(sample_file, T=Float32)
        @test eltype(data_f32.times) == Float32
        @test eltype(data_f32.energies) == Float32

        # Test Int64 conversion
        data_i64 = readevents(sample_file, T=Int64)
        @test eltype(data_i64.times) == Int64
        @test eltype(data_i64.energies) == Int64

        rm(test_dir, recursive=true, force=true)
    end

    @testset "Missing columns" begin
        test_dir = mktempdir()
    
        # Test file with only TIME column
        time_only_file = joinpath(test_dir, "time_only.fits")
        FITS(time_only_file, "w") do f
            write(f, Float32[])
            data = Dict{String,Vector{Float64}}(
                "TIME" => Float64[1.0, 2.0, 3.0]
            )
            write(f, data)
        end
    
        data_time = readevents(time_only_file)
        @test length(data_time.times) == 3
        @test isempty(data_time.energies)
    
        # Test file with only ENERGY column
        energy_only_file = joinpath(test_dir, "energy_only.fits")
        FITS(energy_only_file, "w") do f
            write(f, Float32[])
            data = Dict{String,Vector{Float64}}(
                "ENERGY" => Float64[10.0, 20.0, 30.0]
            )
            write(f, data)
        end
    
        data_energy = readevents(energy_only_file)
        @test isempty(data_energy.times)  # Should be empty as no TIME column exists
        @test length(data_energy.energies) == 3  # ENERGY data should be read even if TIME is missing
    
        rm(test_dir, recursive=true, force=true)
    end
    @testset "Multiple HDUs with TIME columns" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "multiple_hdus.fits")
        
        FITS(sample_file, "w") do f
            write(f, Float32[])
            
            # First HDU with TIME and ENERGY
            data1 = Dict{String,Vector{Float64}}(
                "TIME" => Float64[1:10...],
                "ENERGY" => Float64[11:20...]
            )
            write(f, data1)
            
            # Second HDU with TIME only (should be ignored)
            data2 = Dict{String,Vector{Float64}}(
                "TIME" => Float64[21:30...]
            )
            write(f, data2)
        end

        event_list = readevents(sample_file)
        @test length(event_list.times) == 10
        @test length(event_list.energies) == 10
        @test event_list.times == collect(1:10)
        @test event_list.energies == collect(11:20)
        @test length(event_list.metadata.headers) == 2  # Should only include headers up to first event HDU

        rm(test_dir, recursive=true, force=true)
    end

    @testset "Real data files" begin
        test_filepath = joinpath("data", "monol_testA.evt")
        if isfile(test_filepath)
            data = readevents(test_filepath)
            @test data.filename == test_filepath
            @test length(data.metadata.headers) > 0
            @test !isempty(data.times)
        else
            @info "Test file 'monol_testA.evt' not found. Skipping this test."
        end
    end

    @testset "Error handling" begin
        # Test with non-existent file
        @test_throws Exception readevents("non_existent_file.fits")

        # Test with invalid FITS file
        invalid_file = tempname()
        write(invalid_file, "This is not a FITS file")
        @test_throws Exception readevents(invalid_file)
        rm(invalid_file, force=true)
    end
end