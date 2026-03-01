
@testset "Plotting Recipes" begin
    # Setup: Create a standard EventList and LightCurve for testing
    times = collect(1.0:0.1:10.0)
    energies = rand(length(times)) .* 100
    fake_gti = [2.0 8.0]
    meta = FITSMetadata("test.fits", 1, "keV", Dict{String,Vector}(), Dict{String,Any}(), fake_gti, "GTI", nothing, nothing, nothing, nothing, nothing, nothing)
    el = EventList(times, energies, meta)
    lc = create_lightcurve(el, 1.0)

    @testset "EventList Recipes" begin
        # Test basic plot
        p = plot(el, 0.5)
        @test p isa Plots.Plot
        
        # Test with GTI shading enabled
        p_gti = plot(el, 0.5, show_gtis=true)
        @test p_gti isa Plots.Plot
        
        # Test energy filtering logic
        @test_nowarn plot(el, 0.5, energy_filter=(10.0, 50.0))
    end

    @testset "LightCurve Recipes" begin
        # Test basic LC plot with errors
        p = plot(lc, show_errors=true)
        @test p isa Plots.Plot
        
        # Test robustness: ensure it doesn't crash if GTI is missing
        lc_no_gti = create_lightcurve(EventList(times, energies), 1.0)
        @test_nowarn plot(lc_no_gti)
    end
end
