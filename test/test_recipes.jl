using Test
using Plots

# EventList plotting tests
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0) isa Plots.Plot
    @test plot(el, 2.0) isa Plots.Plot
end

# Empty EventList error handling
let
    empty_el = EventList(Float64[], Float64[])
    @test_throws ErrorException plot(empty_el, 1.0)
end

# Time range filtering
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0, tstart=2.0, tstop=5.0) isa Plots.Plot
    @test plot(el, 1.0, tstart=2.0) isa Plots.Plot
    @test plot(el, 1.0, tstop=4.0) isa Plots.Plot
end

# Energy filtering
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0, energy_filter=(15.0, 45.0)) isa Plots.Plot
    @test plot(el, 1.0, energy_filter=(10.0, 30.0)) isa Plots.Plot
end

# Error display options
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0, show_errors=true) isa Plots.Plot
    @test plot(el, 1.0, show_errors=false) isa Plots.Plot
    @test plot(el, 1.0, show_errors=true, err_method=:poisson) isa Plots.Plot
end

# GTI visualization
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    gti_matrix = [1.0 3.0; 4.0 6.0]
    meta = FITSMetadata{Dict{String,Any}}(
        "test.fits", 2, "ENERGY", Dict{String,Vector}(), Dict{String,Any}(),
        gti_matrix, "GTI"
    )
    el = EventList(times, energies, meta)
    
    @test plot(el, 1.0, show_gtis=true) isa Plots.Plot
    @test plot(el, 1.0, show_gti=true) isa Plots.Plot
    @test plot(el, 1.0, show_gtis=true, gti_alpha=0.5) isa Plots.Plot
end

# BTI visualization
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    gti_matrix = [1.0 3.0; 4.0 6.0]
    meta = FITSMetadata{Dict{String,Any}}(
        "test.fits", 2, "ENERGY", Dict{String,Vector}(), Dict{String,Any}(),
        gti_matrix, "GTI"
    )
    el = EventList(times, energies, meta)
    
    @test plot(el, 1.0, show_btis=true) isa Plots.Plot
    @test plot(el, 1.0, show_bti=true) isa Plots.Plot
    @test plot(el, 1.0, show_btis=true, bti_alpha=0.4) isa Plots.Plot
end

# External GTI support
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    external_gtis = [1.0 2.5; 3.5 5.0]
    
    @test plot(el, 1.0, show_gtis=true, gtis=external_gtis) isa Plots.Plot
    @test plot(el, 1.0, show_btis=true, gtis=external_gtis) isa Plots.Plot
end

# Axis limits configuration
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0, axis_limits=[1.5, 4.5, 0, 10]) isa Plots.Plot
    @test plot(el, 1.0, axis_limits=[1.5, 4.5]) isa Plots.Plot
    @test plot(el, 1.0, axis_limits=[nothing, 4.5, 0, nothing]) isa Plots.Plot
end

# Invalid axis limits warning
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test_logs (:info, r"Created light curve") (:warn, "axis_limits should be a vector of length 2 or 4: [xmin, xmax] or [xmin, xmax, ymin, ymax]") plot(el, 1.0, axis_limits=[1.5, 4.5, 0])
end

# LightCurve plotting
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    errors = [3.0, 4.0, 5.0, 3.5, 2.8]
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc = LightCurve(times, 1.0, counts, errors, nothing, EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc) isa Plots.Plot
    @test plot(lc, show_errors=true) isa Plots.Plot
    @test plot(lc, show_errors=false) isa Plots.Plot
end

# LightCurve without errors
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc = LightCurve(times, 1.0, counts, nothing, nothing, EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc, show_errors=true) isa Plots.Plot
end

# LightCurve with event properties
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    errors = [3.0, 4.0, 5.0, 3.5, 2.8]
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    properties = [EventProperty(:mean_energy, [25.0, 30.0, 35.0, 28.0, 22.0], "keV")]
    lc = LightCurve(times, 1.0, counts, errors, nothing, properties, metadata, :poisson)
    
    @test plot(lc, show_properties=true) isa Plots.Plot
    @test plot(lc, show_properties=true, property_name=:mean_energy) isa Plots.Plot
end

# Non-existent property handling
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    errors = [3.0, 4.0, 5.0, 3.5, 2.8]
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc = LightCurve(times, 1.0, counts, errors, nothing, EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc, show_properties=true, property_name=:nonexistent) isa Plots.Plot
end

# LightCurve with GTI metadata
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    errors = [3.0, 4.0, 5.0, 3.5, 2.8]
    extra_data = Dict{String,Any}("gti_applied" => true, "gti_bounds" => [1.0, 5.0])
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], extra_data)
    lc = LightCurve(times, 1.0, counts, errors, nothing, EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc, show_gtis=true) isa Plots.Plot
    @test plot(lc, show_btis=true) isa Plots.Plot
end

# LightCurve axis limits
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    counts = [10, 15, 20, 12, 8]
    errors = [3.0, 4.0, 5.0, 3.5, 2.8]
    metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 5.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc = LightCurve(times, 1.0, counts, errors, nothing, EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc, axis_limits=[1.5, 4.5, 5, 25]) isa Plots.Plot
    @test plot(lc, axis_limits=[1.5, 4.5]) isa Plots.Plot
end

# Segmented LightCurve plotting
let
    times1 = [1.0, 2.0, 3.0]
    counts1 = [10, 15, 20]
    errors1 = [3.0, 4.0, 5.0]
    metadata1 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 3.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc1 = LightCurve(times1, 1.0, counts1, errors1, nothing, EventProperty{Float64}[], metadata1, :poisson)
    
    times2 = [4.0, 5.0, 6.0]
    counts2 = [12, 8, 18]
    errors2 = [3.5, 2.8, 4.2]
    metadata2 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (4.0, 6.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc2 = LightCurve(times2, 1.0, counts2, errors2, nothing, EventProperty{Float64}[], metadata2, :poisson)
    
    segments = [lc1, lc2]
    
    @test plot(segments) isa Plots.Plot
    @test plot(segments, show_errors=true) isa Plots.Plot
    @test plot(segments, show_segment_boundaries=true) isa Plots.Plot
    @test plot(segments, show_segment_boundaries=false) isa Plots.Plot
end

# Segmented LightCurve with custom colors
let
    times1 = [1.0, 2.0, 3.0]
    counts1 = [10, 15, 20]
    errors1 = [3.0, 4.0, 5.0]
    metadata1 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 3.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc1 = LightCurve(times1, 1.0, counts1, errors1, nothing, EventProperty{Float64}[], metadata1, :poisson)
    
    times2 = [4.0, 5.0, 6.0]
    counts2 = [12, 8, 18]
    errors2 = [3.5, 2.8, 4.2]
    metadata2 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (4.0, 6.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc2 = LightCurve(times2, 1.0, counts2, errors2, nothing, EventProperty{Float64}[], metadata2, :poisson)
    
    segments = [lc1, lc2]
    custom_colors = [:red, :green]
    
    @test plot(segments, segment_colors=custom_colors) isa Plots.Plot
end

# Segmented LightCurve axis limits
let
    times1 = [1.0, 2.0, 3.0]
    counts1 = [10, 15, 20]
    errors1 = [3.0, 4.0, 5.0]
    metadata1 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 3.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc1 = LightCurve(times1, 1.0, counts1, errors1, nothing, EventProperty{Float64}[], metadata1, :poisson)
    
    times2 = [4.0, 5.0, 6.0]
    counts2 = [12, 8, 18]
    errors2 = [3.5, 2.8, 4.2]
    metadata2 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (4.0, 6.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc2 = LightCurve(times2, 1.0, counts2, errors2, nothing, EventProperty{Float64}[], metadata2, :poisson)
    
    segments = [lc1, lc2]
    
    @test plot(segments, axis_limits=[1.5, 5.5, 5, 25]) isa Plots.Plot
    @test plot(segments, axis_limits=[1.5, 5.5]) isa Plots.Plot
end

# Single segment handling
let
    times1 = [1.0, 2.0, 3.0]
    counts1 = [10, 15, 20]
    errors1 = [3.0, 4.0, 5.0]
    metadata1 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 3.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc1 = LightCurve(times1, 1.0, counts1, errors1, nothing, EventProperty{Float64}[], metadata1, :poisson)
    
    segments = [lc1]
    
    @test plot(segments) isa Plots.Plot
    @test plot(segments, show_segment_boundaries=true) isa Plots.Plot
end

# Segmented LightCurve without errors
let
    times1 = [1.0, 2.0, 3.0]
    counts1 = [10, 15, 20]
    metadata1 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (1.0, 3.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc1 = LightCurve(times1, 1.0, counts1, nothing, nothing, EventProperty{Float64}[], metadata1, :poisson)
    
    times2 = [4.0, 5.0, 6.0]
    counts2 = [12, 8, 18]
    metadata2 = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (4.0, 6.0), 1.0, Dict{String,Any}[], Dict{String,Any}())
    lc2 = LightCurve(times2, 1.0, counts2, nothing, nothing, EventProperty{Float64}[], metadata2, :poisson)
    
    segments = [lc1, lc2]
    
    @test plot(segments, show_errors=true) isa Plots.Plot
end

# Color cycling for multiple segments
let
    segments = LightCurve[]
    for i in 1:10
        times = [Float64(i), Float64(i+1)]
        counts = [10, 15]
        errors = [3.0, 4.0]
        metadata = LightCurveMetadata("TEST", "TEST_INST", "TEST_OBJ", 58000.0, (Float64(i), Float64(i+1)), 1.0, Dict{String,Any}[], Dict{String,Any}())
        lc = LightCurve(times, 1.0, counts, errors, nothing, EventProperty{Float64}[], metadata, :poisson)
        push!(segments, lc)
    end
    
    @test plot(segments) isa Plots.Plot
end

# GTI file loading error handling
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0]
    el = EventList(times, energies)
    
    @test plot(el, 1.0, show_gtis=true, gti_file="nonexistent.fits") isa Plots.Plot
end

# Complex GTI/BTI visualization
let
    times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    energies = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    gti_matrix = [1.0 2.5; 3.5 5.0; 6.5 8.0]
    meta = FITSMetadata{Dict{String,Any}}(
        "test.fits", 2, "ENERGY", Dict{String,Vector}(), Dict{String,Any}(),
        gti_matrix, "GTI"
    )
    el = EventList(times, energies, meta)
    
    @test plot(el, 1.0, show_gtis=true, show_btis=true, tstart=0.5, tstop=8.5) isa Plots.Plot
end

# GTI boundary edge cases
let
    times = [2.0, 3.0, 4.0]
    energies = [10.0, 20.0, 30.0]
    gti_matrix = [1.0 2.5; 3.5 5.0]
    meta = FITSMetadata{Dict{String,Any}}(
        "test.fits", 2, "ENERGY", Dict{String,Vector}(), Dict{String,Any}(),
        gti_matrix, "GTI"
    )
    el = EventList(times, energies, meta)
    
    @test plot(el, 1.0, show_gtis=true, tstart=2.0, tstop=4.0) isa Plots.Plot
    @test plot(el, 1.0, show_btis=true, tstart=2.0, tstop=4.0) isa Plots.Plot
end
# Basic rebinning functionality test
let
    times = collect(1.0:0.1:10.0)
    counts = rand(1:10, length(times))
    
    # Create mock light curve
    metadata = LightCurveMetadata(
        "TEST", "TEST", "TEST", 0.0, (1.0, 10.0), 0.1,
        [Dict{String,Any}()], Dict{String,Any}()
    )
    lc = LightCurve(times, 0.1, counts, nothing, nothing, 
                   EventProperty{Float64}[], metadata, :poisson)
    
    @test plot(lc, 1.0) isa Plots.Plot
    @test plot(lc, 0.5) isa Plots.Plot
    @test plot(lc, 2.0) isa Plots.Plot
end