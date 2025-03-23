@testset "LightCurve Tests" begin
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "sample_lightcurve.fits")
    
    # Sample Data
    times = [0.1, 0.5, 1.2, 2.4, 3.0, 4.1, 5.6, 6.8, 7.9, 9.2]
    energies = [100, 200, 150, 180, 250, 300, 275, 220, 190, 210]
    metadata = Dict("TELESCOP" => "TestScope", "MISSION" => "TestMission")
    
    # Create a sample FITS file for testing
    f = FITS(sample_file, "w")
    write(f, Int[])  # Empty primary array
    table = Dict("TIME" => times, "ENERGY" => energies)
    write(f, table)
    close(f)

    @test isfile(sample_file)

    # Test reading the sample file
    eventlist = readevents(sample_file)  # No module prefix since it's directly available
    @test eventlist.filename == sample_file
    @test length(eventlist.times) == length(times)
    @test length(eventlist.energies) == length(energies)

    # Create a LightCurve
    bin_size = 1.0
    err_method = :poisson
    lightcurve = create_lightcurve(eventlist, bin_size, err_method=err_method)

    # Validate LightCurve properties
    @test all(lightcurve.counts .>= 0)
    @test lightcurve.err_method == err_method

    println("LightCurve Tests Passed!")
end