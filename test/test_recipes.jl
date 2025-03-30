#just playing with graphs:) i liked my implmentaion eheheh
@testset "Plotting Recipes" begin
    @testset "EventList Light Curve Recipe" begin
        # Create sample data for testing
        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        energies = [1.0, 2.0, 3.0, 4.0, 5.0]
        events = EventList("test.fits", times, energies, DictMetadata([]))
        
        @testset "Basic Light Curve Plot" begin
            plt = plot(events, 1.0)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Errors" begin
            plt = plot(events, 1.0, show_errors=true)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Gaps" begin
            plt = plot(events, 1.0, show_gaps=true, gap_threshold=2.0)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Axis Limits" begin
            plt = plot(events, 1.0, axis_limits=[0.0, 5.0, 0.0, 3.0])
            @test !isempty(plt.series_list)
        end
    end

    @testset "LightCurve Plot Recipe" begin
        # Create sample data for testing
        timebins = [1.0, 2.0, 3.0, 4.0, 5.0]
        counts = [1, 2, 3, 2, 1]
        count_error = sqrt.(counts)
        lc = LightCurve(timebins, counts, count_error, :poisson)
        
        @testset "Basic Light Curve Plot" begin
            plt = plot(lc)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Errors" begin
            plt = plot(lc, show_errors=true)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Gaps" begin
            plt = plot(lc, show_gaps=true, gap_threshold=2.0)
            @test !isempty(plt.series_list)
        end
        
        @testset "Light Curve Plot with Axis Limits" begin
            plt = plot(lc, axis_limits=[0.0, 5.0, 0.0, 4.0])
            @test !isempty(plt.series_list)
        end
    end

    @testset "EventList Timeline Recipe" begin
        # Create sample data for testing
        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        energies = [1.0, 2.0, 3.0, 4.0, 5.0]
        events = EventList("test.fits", times, energies, DictMetadata([]))
        
        @testset "Basic Event Timeline Plot" begin
            plt = plot(events, Val{:events})
            @test !isempty(plt.series_list)
        end
        
        @testset "Event Timeline Plot with Axis Limits" begin
            plt = plot(events, Val{:events}, axis_limits=[0.0, 5.0, 0.0, 1.0])
            @test !isempty(plt.series_list)
        end
        
        @testset "Event Timeline Plot with Color by Energy" begin
            plt = plot(events, Val{:events}, color_by_energy=true)
            @test !isempty(plt.series_list)
        end
    end

    @testset "Combined LightCurve with Events Recipe" begin
        # Create sample data for testing
        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        energies = [1.0, 2.0, 3.0, 4.0, 5.0]
        events = EventList("test.fits", times, energies, DictMetadata([]))
        timebins = [1.0, 2.0, 3.0, 4.0, 5.0]
        counts = [1, 2, 3, 2, 1]
        count_error = sqrt.(counts)
        lc = LightCurve(timebins, counts, count_error, :poisson)
        
        @testset "Basic Combined Plot" begin
            plt = plot(Val{:lightcurve_with_events}, lc, events)
            @test !isempty(plt.series_list)
        end
        
        @testset "Combined Plot with Axis Limits" begin
            plt = plot(Val{:lightcurve_with_events}, lc, events, axis_limits=[0.0, 5.0, 0.0, 4.0])
            @test !isempty(plt.series_list)
        end
        
        @testset "Combined Plot with Color by Energy" begin
            plt = plot(Val{:lightcurve_with_events}, lc, events, color_by_energy=true)
            @test !isempty(plt.series_list)
        end
    end
end
#this are some basic tests for the plotting recipes in the Stingray.jl package.
#they are not exhaustive and only test the basic functionality of the recipes.
#more tests should be added to cover all the edge cases and possible inputs.
#the tests are written in a way that they can be easily extended and modified.