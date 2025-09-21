# BTI Histogram Recipe Tests
let
    # Suppress prints output during tests(prints to avoid length limit)
    original_stdout = stdout
    redirect_stdout(devnull)
    
    try
        # Test 1: No BTIs
        times = collect(0.0:1.0:100.0)
        energies = rand(length(times)) .* 100
        gti_full = reshape([-1.0, 101.0], 1, 2)
        metadata = FITSMetadata("test_file.fits", 1, "keV", Dict{String,Vector}(), 
                               Dict{String,Any}(), gti_full, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_no_btis = EventList(times, energies, metadata)
        p1 = plot(BTIAnalysisPlot(el_no_btis), bti_analysis=true)
        @test p1 isa Plots.Plot
        
        # Test 2: With BTIs
        times2 = collect(0.0:0.5:50.0)
        energies2 = rand(length(times2)) .* 100
        gti_gaps = [0.0 10.0; 20.0 30.0; 40.0 50.0]  # Creates BTIs at [10-20], [30-40]
        metadata2 = FITSMetadata("test_file2.fits", 1, "keV", Dict{String,Vector}(), 
                                Dict{String,Any}(), gti_gaps, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_with_btis = EventList(times2, energies2, metadata2)
        p2 = plot(BTIAnalysisPlot(el_with_btis), bti_analysis=true)
        @test p2 isa Plots.Plot
        
        # Test 3: Parameter variations with BTIs
        p3 = plot(BTIAnalysisPlot(el_with_btis), bti_analysis=true, nbins=50)
        @test p3 isa Plots.Plot
        p4 = plot(BTIAnalysisPlot(el_with_btis), bti_analysis=true, nbins=20, xlims_range=(5, 15))
        @test p4 isa Plots.Plot
        p5 = plot(BTIAnalysisPlot(el_with_btis), bti_analysis=true, ylims_range=(0.5, 5))
        @test p5 isa Plots.Plot
        
        # Test 4: Complex GTI pattern with multiple BTIs
        times3 = collect(0.0:0.1:30.0)
        energies3 = rand(length(times3)) .* 100
        # Creates multiple BTIs: [5-6], [10-12], [15-18], [22-25]
        gti_complex = [0.0 5.0; 6.0 10.0; 12.0 15.0; 18.0 22.0; 25.0 30.0]
        metadata3 = FITSMetadata("test_file3.fits", 1, "keV", Dict{String,Vector}(),
                                Dict{String,Any}(), gti_complex, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_complex = EventList(times3, energies3, metadata3)
        p6 = plot(BTIAnalysisPlot(el_complex), bti_analysis=true)
        @test p6 isa Plots.Plot
        
        # Test 5: Single very large BTI
        times4 = [0.0, 1.0, 100.0, 101.0]
        energies4 = [10.0, 20.0, 30.0, 40.0]
        gti_large_gap = [0.0 1.0; 100.0 101.0]  # Creates 99s BTI gap
        metadata4 = FITSMetadata("test_file4.fits", 1, "keV", Dict{String,Vector}(),
                                Dict{String,Any}(), gti_large_gap, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_large_gap = EventList(times4, energies4, metadata4)
        p8 = plot(BTIAnalysisPlot(el_large_gap), bti_analysis=true)
        @test p8 isa Plots.Plot
        
        # Test 6: Edge case - very short BTIs only
        times5 = collect(0.0:0.01:10.0)
        energies5 = rand(length(times5)) .* 100
        # Creates very short BTIs: [1-1.1], [2-2.05], [3-3.02]
        gti_short = [0.0 1.0; 1.1 2.0; 2.05 3.0; 3.02 10.0]
        metadata5 = FITSMetadata("test_file5.fits", 1, "keV", Dict{String,Vector}(),
                                Dict{String,Any}(), gti_short, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_short_btis = EventList(times5, energies5, metadata5)
        p9 = plot(BTIAnalysisPlot(el_short_btis), bti_analysis=true)
        @test p9 isa Plots.Plot
        
        # Test 7: Mixed BTI lengths (short and long)
        times6 = collect(0.0:0.1:100.0)
        energies6 = rand(length(times6)) .* 100
        # Creates mixed BTIs: [5-5.1] (0.1s), [10-20] (10s), [30-30.5] (0.5s), [50-80] (30s)
        gti_mixed = [0.0 5.0; 5.1 10.0; 20.0 30.0; 30.5 50.0; 80.0 100.0]
        metadata6 = FITSMetadata("test_file6.fits", 1, "keV", Dict{String,Vector}(),
                                Dict{String,Any}(), gti_mixed, "GTI",nothing, nothing, nothing, nothing, nothing, nothing)
        el_mixed_btis = EventList(times6, energies6, metadata6)
        p10 = plot(BTIAnalysisPlot(el_mixed_btis), bti_analysis=true)
        @test p10 isa Plots.Plot
        
    finally
        # Always restore stdout
        redirect_stdout(original_stdout)
    end
end