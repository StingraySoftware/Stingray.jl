@testset "LightCurve Implementation Tests" begin
    @testset "Structure Tests" begin
        # Test EventProperty structure
        @testset "EventProperty" begin
            prop = EventProperty{Float64}(:test, [1.0, 2.0, 3.0], "units")
            @test prop.name === :test
            @test prop.values == [1.0, 2.0, 3.0]
            @test prop.unit == "units"
            @test typeof(prop) <: EventProperty{Float64}
        end

        # Test LightCurveMetadata structure
        @testset "LightCurveMetadata" begin
            metadata = LightCurveMetadata(
                "TEST_TELESCOPE",
                "TEST_INSTRUMENT",
                "TEST_OBJECT",
                58000.0,
                (0.0, 100.0),
                1.0,
                [Dict{String,Any}("TEST" => "VALUE")],
                Dict{String,Any}("extra_info" => "test")
            )
            @test metadata.telescope == "TEST_TELESCOPE"
            @test metadata.instrument == "TEST_INSTRUMENT"
            @test metadata.object == "TEST_OBJECT"
            @test metadata.mjdref == 58000.0
            @test metadata.time_range == (0.0, 100.0)
            @test metadata.bin_size == 1.0
            @test length(metadata.headers) == 1
            @test haskey(metadata.extra, "extra_info")
            @test metadata.extra["extra_info"] == "test"
        end

        # Test LightCurve structure
        @testset "LightCurve Basic Structure" begin
            timebins = [1.5, 2.5, 3.5]
            bin_edges = [1.0, 2.0, 3.0, 4.0]
            counts = [1, 2, 1]
            errors = Float64[1.0, √2, 1.0]
            exposure = fill(1.0, 3)
            props = [EventProperty{Float64}(:test, [1.0, 2.0, 3.0], "units")]
            metadata = LightCurveMetadata(
                "TEST", "TEST", "TEST", 0.0, (1.0, 4.0), 1.0,
                [Dict{String,Any}()], Dict{String,Any}()
            )
            
            lc = LightCurve{Float64}(
                timebins, bin_edges, counts, errors, exposure,
                props, metadata, :poisson
            )
            
            @test lc.timebins == timebins
            @test lc.bin_edges == bin_edges
            @test lc.counts == counts
            @test lc.count_error == errors
            @test lc.exposure == exposure
            @test length(lc.properties) == 1
            @test lc.err_method === :poisson
            @test typeof(lc) <: AbstractLightCurve{Float64}
        end
    end

    @testset "Error Calculation Tests" begin
        @testset "Error Methods" begin
            # Test Poisson errors
            counts = [0, 1, 4, 9, 16]
            exposure = fill(1.0, length(counts))
            
            errors = calculate_errors(counts, :poisson, exposure)
            @test errors ≈ [1.0, 1.0, 2.0, 3.0, 4.0]
            
            # Test Gaussian errors
            gaussian_errs = [0.5, 1.0, 1.5, 2.0, 2.5]
            errors_gauss = calculate_errors(counts, :gaussian, exposure, 
                                         gaussian_errors=gaussian_errs)
            @test errors_gauss == gaussian_errs
            
            # Test error conditions
            @test_throws ArgumentError calculate_errors(counts, :gaussian, exposure)
            @test_throws ArgumentError calculate_errors(
                counts, :gaussian, exposure, 
                gaussian_errors=[1.0, 2.0]
            )
            @test_throws ArgumentError calculate_errors(counts, :invalid, exposure)
        end
    end

    @testset "Input Validation" begin
        @testset "validate_lightcurve_inputs" begin
            # Test valid inputs
            valid_events = EventList{Float64}(
                "test.fits",
                [1.0, 2.0, 3.0],
                [10.0, 20.0, 30.0],
                Dict{String,Vector}(),
                DictMetadata([Dict{String,Any}()])
            )
            
            @test_nowarn validate_lightcurve_inputs(valid_events, 1.0, :poisson, nothing)
            
            # Test invalid bin size
            @test_throws ArgumentError validate_lightcurve_inputs(valid_events, 0.0, :poisson, nothing)
            @test_throws ArgumentError validate_lightcurve_inputs(valid_events, -1.0, :poisson, nothing)
            
            # Test invalid error method
            @test_throws ArgumentError validate_lightcurve_inputs(valid_events, 1.0, :invalid, nothing)
            
            # Test missing gaussian errors
            @test_throws ArgumentError validate_lightcurve_inputs(valid_events, 1.0, :gaussian, nothing)
        end

        @testset "Event Filtering" begin
            times = [1.0, 2.0, 3.0, 4.0, 5.0]
            energies = [10.0, 20.0, 30.0, 40.0, 50.0]
            
            # Test time filtering
            filtered_times, filtered_energies, start_t, stop_t = 
                apply_event_filters(times, energies, 2.0, 4.0, nothing)
            @test all(2.0 .<= filtered_times .<= 4.0)
            @test length(filtered_times) == 3
            @test start_t == 2.0
            @test stop_t == 4.0
            
            # Test energy filtering
            filtered_times, filtered_energies, start_t, stop_t = 
                apply_event_filters(times, energies, nothing, nothing, (15.0, 35.0))
            @test all(15.0 .<= filtered_energies .< 35.0)
            
            # Test combined filtering
            filtered_times, filtered_energies, start_t, stop_t = 
                apply_event_filters(times, energies, 2.0, 4.0, (15.0, 35.0))
            @test all(2.0 .<= filtered_times .<= 4.0)
            @test all(15.0 .<= filtered_energies .< 35.0)
        end
    end

    @testset "Binning Operations" begin
        @testset "Time Bin Creation" begin
            start_time = 1.0
            stop_time = 5.0
            binsize = 1.0
            
            edges, centers = create_time_bins(start_time, stop_time, binsize)
            num_bins = ceil(Int, (stop_time - start_time) / binsize)
            
            expected_edges = [start_time + i * binsize for i in 0:(num_bins)]
            expected_centers = [start_time + (i + 0.5) * binsize for i in 0:(num_bins-1)]
            
            @test length(edges) == length(expected_edges)
            @test length(centers) == length(expected_centers)
            @test all(isapprox.(edges, expected_edges, rtol=1e-10))
            @test all(isapprox.(centers, expected_centers, rtol=1e-10))
            
            # Test with fractional boundaries
            edges_frac, centers_frac = create_time_bins(0.5, 2.5, 0.5)
            @test isapprox(edges_frac[1], 0.5, rtol=1e-10)
            @test edges_frac[end] >= 2.5
            @test isapprox(centers_frac[1], 0.75, rtol=1e-10)
        end

        @testset "Event Binning" begin
            times = [1.1, 1.2, 2.3, 2.4, 3.5]
            edges = [1.0, 2.0, 3.0, 4.0]
            
            counts = bin_events(times, edges)
            @test counts == [2, 2, 1]
            
            # Test empty data
            @test all(bin_events(Float64[], edges) .== 0)
            
            # Test single event
            @test bin_events([1.5], edges) == [1, 0, 0]
        end
    end

    @testset "Property Calculations" begin
        @testset "Additional Properties" begin
            times = [1.1, 1.2, 2.3, 2.4, 3.5]
            energies = [10.0, 20.0, 15.0, 25.0, 30.0]
            edges = [1.0, 2.0, 3.0, 4.0]
            centers = [1.5, 2.5, 3.5]
            
            props = calculate_additional_properties(
                times, energies, edges, centers
            )
            
            @test length(props) == 1
            @test props[1].name === :mean_energy
            @test props[1].unit == "keV"
            @test length(props[1].values) == length(centers)
            
            # Test mean energy calculation
            mean_energies = props[1].values
            @test mean_energies[1] ≈ mean([10.0, 20.0])
            @test mean_energies[2] ≈ mean([15.0, 25.0])
            @test mean_energies[3] ≈ 30.0
            
            # Test without energies
            props_no_energy = calculate_additional_properties(
                times, nothing, edges, centers
            )
            @test isempty(props_no_energy)
        end
    end

    @testset "Rebinning" begin
        @testset "Basic Rebinning" begin
            start_time = 1.0
            end_time = 7.0
            old_binsize = 0.5
            new_binsize = 1.0
            
            # Create times and edges that align perfectly with both bin sizes
            times = collect(start_time + old_binsize/2 : old_binsize : end_time - old_binsize/2)
            edges = collect(start_time : old_binsize : end_time)
            counts = ones(Int, length(times))
            
            lc = LightCurve{Float64}(
                times,
                edges,
                counts,
                sqrt.(Float64.(counts)),
                fill(old_binsize, length(times)),
                Vector{EventProperty{Float64}}(),
                LightCurveMetadata(
                    "TEST", "TEST", "TEST", 0.0,
                    (start_time, end_time), old_binsize,
                    [Dict{String,Any}()],
                    Dict{String,Any}()
                ),
                :poisson
            )
            
            # Test rebinning to larger bins
            new_lc = rebin(lc, new_binsize)
            
            # Calculate expected number of bins
            expected_bins = ceil(Int, (end_time - start_time) / new_binsize)
            @test length(new_lc.counts) == expected_bins
            @test all(new_lc.exposure .== new_binsize)
            @test sum(new_lc.counts) == sum(lc.counts)
        end

        @testset "Property Rebinning" begin
            start_time = 1.0
            end_time = 7.0
            old_binsize = 1.0
            new_binsize = 2.0
            
            times = collect(start_time + old_binsize/2 : old_binsize : end_time - old_binsize/2)
            edges = collect(start_time : old_binsize : end_time)
            n_bins = length(times)
            
            counts = fill(2, n_bins)
            energy_values = collect(10.0:10.0:(10.0*n_bins))
            props = [EventProperty{Float64}(:mean_energy, energy_values, "keV")]
            
            lc = LightCurve{Float64}(
                times,
                edges,
                counts,
                sqrt.(Float64.(counts)),
                fill(old_binsize, n_bins),
                props,
                LightCurveMetadata(
                    "TEST", "TEST", "TEST", 0.0,
                    (start_time, end_time), old_binsize,
                    [Dict{String,Any}()],
                    Dict{String,Any}()
                ),
                :poisson
            )
            
            # Test rebinning with exact factor
            new_lc = rebin(lc, new_binsize)
            
            start_bin = floor(start_time / new_binsize) * new_binsize
            num_new_bins = ceil(Int, (end_time - start_bin) / new_binsize)
            
            @test new_lc.metadata.bin_size == new_binsize
            @test sum(new_lc.counts) == sum(lc.counts)
            @test length(new_lc.properties) == length(lc.properties)
            @test all(new_lc.exposure .== new_binsize)
            
            # Test half range rebinning
            total_range = end_time - start_time
            half_range_size = total_range / 2
            lc_half = rebin(lc, half_range_size)
            
            start_half = floor(start_time / half_range_size) * half_range_size
            n_half_bins = ceil(Int, (end_time - start_half) / half_range_size)
            @test length(lc_half.counts) == n_half_bins
            @test sum(lc_half.counts) == sum(lc.counts)
        end
    end

    @testset "Array Interface" begin
        times = [1.5, 2.5, 3.5]
        counts = [1, 2, 1]
        lc = LightCurve{Float64}(
            times,
            [1.0, 2.0, 3.0, 4.0],
            counts,
            sqrt.(Float64.(counts)),
            fill(1.0, 3),
            Vector{EventProperty{Float64}}(),
            LightCurveMetadata(
                "TEST", "TEST", "TEST", 0.0,
                (1.0, 4.0), 1.0,
                [Dict{String,Any}()],
                Dict{String,Any}()
            ),
            :poisson
        )
        
        @test length(lc) == 3
        @test size(lc) == (3,)
        @test lc[1] == (1.5, 1)
        @test lc[2] == (2.5, 2)
        @test lc[3] == (3.5, 1)
    end
end
