@testset "LightCurve Tests" begin
    # Test 1: Create a lightcurve from a sample EventList
    @testset "Lightcurve creation" begin
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
        
        # Create a lightcurve from the event data
        lc =create_lightcurve(data, 1.0)
        
        # Basic size checks
        @test length(lc.timebins) == 6  # 5 bins + 1 endpoint
        @test length(lc.counts) == 5
        
        # Time bin checks
        @test lc.timebins[1] ≈ 1.0
        @test lc.timebins[end] ≈ 6.0
        
        # Count checks - each bin should have 1 count
        @test all(lc.counts .== 1)
        
        # Error checks - since each bin has 1 count, errors should be sqrt(1)
        @test all(lc.count_error .≈ fill(sqrt(1.0), 5))
        
        # Method check
        @test lc.err_method == :poisson
    end

   # Test 2: Test lightcurve with different bin sizes
    @testset "Different bin sizes" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_bin_sizes.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])
        # Create data
        times = Float64[1.0, 2.0, 3.0, 3.5, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 30.0, 25.0, 15.0, 10.0]
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        data = readevents(sample_file)
    
        # Test with bin size of 2.0
        lc_2 =create_lightcurve(data, 2.0)
        @test length(lc_2.timebins) == 4  # 3 bins + 1 endpoint
        @test length(lc_2.counts) == 3
        @test lc_2.timebins[1] == 1.0
        @test lc_2.timebins[end] == 7.0  # Updated expected value
        @test lc_2.counts == [2, 3, 1]  # Updated expected value
        @test lc_2.count_error == sqrt.([2, 3, 1])  # Updated expected value
    
        # Test with bin size of 1.0
        lc_1 =create_lightcurve(data, 1.0)
        @test length(lc_1.timebins) == 6  # 5 bins + 1 endpoint
        @test length(lc_1.counts) == 5
        @test lc_1.timebins[1] == 1.0
        @test lc_1.timebins[end] == 6.0
        @test lc_1.counts == [1, 1, 2, 1, 1]  # Updated expected value
        @test lc_1.count_error == sqrt.([1, 1, 2, 1, 1])  # Updated expected value
    end
    
    # Test 3: Test lightcurve error computation methods
    @testset "Error computation methods" begin
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_error_methods.fits")
        f = FITS(sample_file, "w")
        write(f, Int[])
        # Create data
        times = Float64[1.0, 2.0, 2.5, 3.0, 4.5, 5.0]
        energies = Float64[10.0, 20.0, 20.0, 25.0, 15.0, 10.0]
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)
        data = readevents(sample_file)
    
        # Test with default error method (Poisson)
        lc =create_lightcurve(data, 1.0)
        @test lc.count_error == sqrt.([1, 2, 1, 1, 1])  # Updated expected value
        
        # Add more error computation methods here if needed
    end

    # Test 4: Struct Type Validation
    @testset "LightCurve Struct Type Checks" begin
        # Create a sample FITS file for type testing
        test_dir = mktempdir()
        sample_file = joinpath(test_dir, "sample_types.fits")
        
        # Prepare test data
        f = FITS(sample_file, "w")
        write(f, Int[])  # Empty primary array
        
        # Create test data
        times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
        energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
        
        table = Dict{String,Array}()
        table["TIME"] = times
        table["ENERGY"] = energies
        write(f, table)
        close(f)

        # Test type-specific instantiations
        @testset "Type Parametric Struct Tests" begin
            # Test Float64 EventList
            data_f64 = readevents(sample_file, T = Float64)
            @test isa(data_f64, EventList{Float64})
            @test typeof(data_f64) == EventList{Float64}
            
            # Test lightcurve creation
            lc_f64 =create_lightcurve(data_f64, 1.0)
            @test isa(lc_f64, LightCurve{Float64})
            @test typeof(lc_f64) == LightCurve{Float64}

            # Test Float32 EventList
            data_f32 = readevents(sample_file, T = Float32)
            @test isa(data_f32, EventList{Float32})
            @test typeof(data_f32) == EventList{Float32}

            # Test lightcurve creation
            lc_f32 =create_lightcurve(data_f32, 1.0)
            @test isa(lc_f32, LightCurve{Float32})
            @test typeof(lc_f32) == LightCurve{Float32}

            # Test Int64 EventList
            data_i64 = readevents(sample_file, T = Int64)
            @test isa(data_i64, EventList{Int64})
            @test typeof(data_i64) == EventList{Int64}

            # Test lightcurve creation
            lc_i64 =create_lightcurve(data_i64, 1.0)
            @test isa(lc_i64, LightCurve{Int64})
            @test typeof(lc_i64) == LightCurve{Int64}
        end

        # Test struct field types
        @testset "Struct Field Type Checks" begin
            data = readevents(sample_file)
            
            # Check filename type
            @test isa(data.filename, String)
            
            # Check times and energies vector types
            @test isa(data.times, Vector{Float64})
            @test isa(data.energies, Vector{Float64})
            
            # Check metadata type
            @test isa(data.metadata, DictMetadata)
            @test isa(data.metadata.headers, Vector{Dict{String,Any}})

            # Test lightcurve creation
            lc =create_lightcurve(data, 1.0)

            # Check timebins type
            @test isa(lc.timebins, Vector{Float64})

            # Check counts type
            @test isa(lc.counts, Vector{Int})

            # Check count_error type
            @test isa(lc.count_error, Vector{Float64})

            # Check err_method type
            @test isa(lc.err_method, Symbol)
        end
    end
end