
# Helper function to create mock FITSMetadata
function create_test_eventlist(times::Vector{Float64}, energies::Union{Vector{Float64}, Nothing}=nothing)
    mock_headers = Dict{String,Any}()
    mock_metadata = FITSMetadata(
        "",  # telescope
        1,   # channel
        nothing, # unit
        Dict{String,Vector}(),  # column info
        mock_headers,  # headers
        [0.0 maximum(times)],  # GTI
        nothing  # maskfile
    )
    
    return EventList(times, energies, mock_metadata)
end

# Helper function to create mock LightCurve with GTI
function create_test_lightcurve(times::Vector{Float64}, counts::Vector{Int}, dt::Float64=1.0)
    metadata = LightCurveMetadata(
        "", "", "", 0.0, 
        (minimum(times)-dt/2, maximum(times)+dt/2), 
        dt, 
        Vector{Dict{String,Any}}(),
        Dict{String,Any}("gti" => [minimum(times) maximum(times)])
    )
    
    return LightCurve(
        times, dt, counts, nothing, fill(dt, length(times)), EventProperty{Float64}[], 
        metadata, :poisson
    )
end


    
let
    times = collect(0.0:0.01:10.0)
    counts = rand(Poisson(100), length(times))
    dt = times[2] - times[1]
    lc = create_test_lightcurve(times, counts, dt)
    # Constructor Tests
    let
        # Test basic construction
        ps = Powerspectrum(lc)
        @test ps !== nothing
        
        # Test with events list
        events = create_test_eventlist(sort(rand(1000) * 10.0))
        lc_from_events = create_lightcurve(events, dt)
        ps_events = Powerspectrum(lc_from_events)
        @test ps_events !== nothing
    end
    # Normalization Tests
    let
        ps = Powerspectrum(lc, norm="leahy")
        @test mean(ps.power) ≈ 2.0 rtol=0.5
        
        ps_frac = Powerspectrum(lc, norm="frac")
        @test all(ps_frac.power .≥ 0)
        
        ps_abs = Powerspectrum(lc, norm="abs")
        @test ps_abs.norm == "abs"
    end
end
# Frequency Properties
let
    # Use exact powers of 2 for length to avoid frequency rounding issues
    n_points = 1024
    dt = 0.01
    times = collect(0.0:dt:(n_points-1)*dt)
    counts = rand(Poisson(100), length(times))
    lc = create_test_lightcurve(times, counts, dt)
    ps = Powerspectrum(lc)
    
    @test issorted(ps.freq)
    @test ps.freq[1] > 0
    @test ps.freq[end] ≤ 1/(2*dt)  # Nyquist frequency
    
    # Test frequency resolution
    tseg = times[end] - times[1] + dt  # Include the last bin width
    expected_df = 1/tseg
    @test abs(ps.df - expected_df) < 1e-10
end
# Error Handling
let
    times = collect(0.0:0.01:10.0)
    counts = rand(Poisson(100), length(times))
    dt = times[2] - times[1]
    lc = create_test_lightcurve(times, counts, dt)
    
    @test_throws ArgumentError AveragedPowerspectrum(lc, 0.0)
    @test_throws ArgumentError AveragedPowerspectrum(lc, -1.0)
    @test_throws ArgumentError Powerspectrum(lc, norm="invalid")
end
#Light Curve Rebinning
let
    times = collect(0.0:0.01:10.0)
    counts = rand(Poisson(100), length(times))
    dt = times[2] - times[1]
    lc = create_test_lightcurve(times, counts, dt)
    
    # Test rebinning with integer factor
    new_dt = 2 * dt
    rebinned_lc = rebin(lc, new_dt)
    
    @test rebinned_lc.metadata.bin_size ≈ new_dt
    @test length(rebinned_lc.time) < length(lc.time)
    @test sum(rebinned_lc.counts) ≤ sum(lc.counts)  # Some counts might be dropped at edges
    
    # Test error handling
    @test_throws ArgumentError rebin(lc, dt/2)  # Cannot decrease bin size
end
# Event List Properties
let
    events = create_test_eventlist(sort(rand(1000) * 10.0))
    dt = 0.1
    
    # Convert to light curve using your create_lightcurve function
    lc_from_events = create_lightcurve(events, dt)
    ps = Powerspectrum(lc_from_events)
    
    @test ps !== nothing
    @test all(ps.freq .≥ 0)
    @test issorted(ps.freq)
    
    # Test averaged power spectrum
    aps = AveragedPowerspectrum(lc_from_events, 1.0)
    @test aps.m > 0
    @test all(aps.power .≥ 0)
end