@testset "Power Spectrum Analysis" begin
    @testset "Basic Spectrum Creation" begin
        dt = 0.001
        T = 100.0
        t = 0:dt:T
        freq = 2.0
        
        signal = 2.0 * sin.(2π * freq * t)
        
        metadata = DictMetadata([Dict{String,Any}(
            "TELESCOP" => "TEST",
            "INSTRUME" => "SYNTHETIC",
            "DATATYPE" => "Single frequency analysis",
        )])
        
        events = EventList{Float64}(
            "test_basic.fits",
            collect(Float64, t),
            Float64.(signal),
            metadata
        )
    
        ps = powerspectrum(events, dt=dt)
        @test isa(ps, AveragedPowerspectrum)  # This is what the function actually returns
        @test length(ps.freqs) == length(ps.power)
        @test length(ps.power) == length(ps.power_errors)
        
        freqs = ps.freqs
        powers = ps.power
        peak_idx = argmax(powers)
        detected_freq = ps.freqs[peak_idx]
        
        @test isapprox(detected_freq, freq, rtol=0.1, atol=1/(T*dt))
    end

    @testset "Multiple Frequency Detection" begin
        # Use power of 2 for better FFT performance
        n_points = 2^17  # Longer time series
        dt = 0.001
        t = range(0, step=dt, length=n_points)
        T = t[end]
        f1, f2 = 1.0, 2.5
        
        # Generate clean signal with high SNR
        signal = 50.0 * (sin.(2π * f1 * t) + sin.(2π * f2 * t))
        
        metadata = DictMetadata([Dict{String,Any}(
            "TELESCOP" => "TEST",
            "INSTRUME" => "SYNTHETIC",
            "DATATYPE" => "Multi-frequency analysis",
        )])
        
        events = EventList(
            "test_multi.fits",
            collect(Float64, t),
            Float64.(signal),
            metadata
        )

        ps = powerspectrum(events, dt=dt)
        
        # Skip DC and very low frequencies
        min_freq_idx = max(2, floor(Int, 0.1/(T*dt)))
        freqs = ps.freqs[min_freq_idx:end]
        powers = ps.power[min_freq_idx:end]
        
        # Find peaks using local maxima and power threshold
        peak_indices = Int[]
        for i in 2:(length(powers)-1)
            if powers[i] > powers[i-1] && powers[i] > powers[i+1]
                push!(peak_indices, i)
            end
        end
        
        if !isempty(peak_indices)
            peak_freqs = freqs[peak_indices]
            peak_powers = powers[peak_indices]
            
            sorted_idx = sortperm(peak_powers, rev=true)
            peak_freqs = peak_freqs[sorted_idx]
            
            df = 1/T 
            found_f1 = any(f -> abs(f - f1) ≤ 8df, peak_freqs)
            found_f2 = any(f -> abs(f - f2) ≤ 8df, peak_freqs)
            
            @info "Peak Detection Results" begin
                expected_freqs = [f1, f2]
                detected_freqs = peak_freqs[1:min(5,end)]
                frequency_resolution = df
                tolerance = 8df
            end
            
            @test found_f1
            @test found_f2
        else
            @test false
        end
    end

     @testset "Averaged Power Spectrum" begin
        n_points = 2^16  # 65536 points
        dt = 0.001
        t = range(0, step=dt, length=n_points)
        freq = 2.0
        
        signal = 20.0 * sin.(2π * freq * t)
        
        metadata = DictMetadata([Dict{String,Any}(
            "TELESCOP" => "TEST",
            "INSTRUME" => "SYNTHETIC",
            "DATATYPE" => "Averaged spectrum analysis"
        )])
        
        events = EventList(
            "test_avg.fits",
            collect(Float64, t),
            Float64.(signal),
            metadata
        )
    
        # Use smaller segments for better averaging
        segment_size = 2048
        
        # Calculate averaged power spectrum
        aps = powerspectrum(events, dt=dt, segment_size=segment_size)
        
        # Basic structure tests
        @test isa(aps, AveragedPowerspectrum)
        @test length(aps.freqs) == length(aps.power)
        @test length(aps.power) == length(aps.power_errors)
        @test aps.m > 0
        
        # Analyze the spectrum
        powers = aps.power
        freqs = aps.freqs
        
        # Find peaks
        peak_idx = argmax(powers)
        detected_freq = freqs[peak_idx]
        
        # Frequency resolution
        df = 1/(segment_size * dt)
        
        # Test with appropriate tolerance
        @test abs(detected_freq - freq) <= 2*df  # Allow 2 frequency bins of tolerance
    end
    @testset "Normalization Methods" begin
        n_points = 2^12
        dt = 0.001
        t = range(0, step=dt, length=n_points)
        signal = sin.(2π * 2.0 * t)
        
        metadata = DictMetadata([Dict{String,Any}(
            "TELESCOP" => "TEST",
            "INSTRUME" => "SYNTHETIC",
            "DATATYPE" => "Normalization analysis",
        )])
        
        events = EventList(
            "test_norm.fits",
            collect(Float64, t),
            Float64.(signal),
            metadata
        )

        for norm in ["leahy", "frac", "abs"]
            ps = powerspectrum(events, dt=dt, norm=norm)
            @test all(ps.power .>= 0)
            @test all(ps.power_errors .>= 0)
            @test length(ps.freqs) == length(ps.power)
        end
    end
end