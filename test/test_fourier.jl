function compare_tables(table1, table2; rtol=0.001, discard = [])

    s_discard = Symbol.(discard)
    test_result = true

    for key in propertynames(table1)
        if key in s_discard
            continue
        end
        oe, oc = getproperty(table1,key), getproperty(table1,key)
        if oe isa Integer || oe isa String
            if !(oe==oc) 
                test_result = false
                break
            end
        elseif isnothing(oe)
            if !(isnothing(oc)) 
                test_result = false
                break
            end
        else
            if !(≈(oe,oc,rtol=rtol)) 
                test_result = false
                break
            end
        end
    end

    table1 = Metadata.drop_metadata(table1)
    table2 = Metadata.drop_metadata(table2)

    for field in names(table1)
        if field in discard
            continue
        end
        oe, oc = table1[!,field], table2[!,field]

        if !(≈(oe,oc,rtol=rtol)) 
            test_result = false
            break
        end
    end
    @test test_result
end

@testset "positive_fft_bins" begin
    freq = fftfreq(11)
    goodbins = positive_fft_bins(11)
    @test filter(x -> x >0, freq) == freq[goodbins]
end

@testset "test_coherence" begin

    fid = h5open(joinpath(@__DIR__ ,"data","sample_variable_lc.h5"),"r")
    data = (fid["dataset1"][1:10000]) .* 1000.0
    close(fid)

    data1 = rand.(Poisson.(data))
    data2 = rand.(Poisson.(data))
    ft1 = fft(data1)
    ft2 = fft(data2)
    dt = 0.01
    N = length(data)
    mean = Statistics.mean(data)
    meanrate = mean / dt
    freq = fftfreq(length(data), 1/dt)
    good = 0 .< freq .< 0.1
    keepat!(ft1, good)
    keepat!(ft2, good)
    cross = normalize_periodograms(
        ft1 .* conj(ft2), dt, N, mean_flux = mean, norm="abs", power_type="all")
    pds1 = normalize_periodograms(
        ft1 .* conj(ft1), dt, N, mean_flux = mean, norm="abs", power_type="real")
    pds2 = normalize_periodograms(
        ft2 .* conj(ft2), dt, N, mean_flux = mean, norm="abs", power_type="real")

    p1noise = poisson_level("abs",meanrate=meanrate)
    p2noise = poisson_level("abs",meanrate=meanrate)

    @testset "test_intrinsic_coherence" begin
        coh = estimate_intrinsic_coherence.(
            cross, pds1, pds2, p1noise, p2noise, N)
        @test all(x -> isapprox(x, 1; atol=0.001), coh)
    end

    @testset "test_raw_high_coherence" begin
        coh = raw_coherence.(cross, pds1, pds2, p1noise, p2noise, N)
        @test all(x -> isapprox(x, 1; atol=0.001), coh)
    end

    @testset "test_raw_low_coherence" begin
        nbins = 2
        C, P1, P2 = @view(cross[1:nbins]), @view(pds1[1:nbins]), @view(pds2[1:nbins])
        bsq = bias_term.(P1, P2, p1noise, p2noise, N)
        # must be lower than bsq!
        low_coh_cross = @. rand(Normal(bsq^0.5 / 10, bsq^0.5 / 100)) + 0.0im
        coh = raw_coherence.(low_coh_cross, P1, P2, p1noise, p2noise, N)
        @test all(iszero,coh)
    end

    @testset "test_raw_high_bias" begin
        C = [12986.0 + 8694.0im]
        P1 = [476156.0]
        P2 = [482751.0]
        P1noise = 495955
        P2noise = 494967
    
        coh = raw_coherence.(C, P1, P2, P1noise, P2noise, 499; intrinsic_coherence=1)
        coh_sngl = raw_coherence(C[1], P1[1], P2[1], P1noise, P2noise, 499; intrinsic_coherence=1)
        @test coh==real(C .* conj(C))./ (P1 .* P2)
        @test coh_sngl==real(C .* conj(C))[1] / (P1[1] * P2[1])
    end
    
end

@testset "test_fourier" begin
    dt = 1
    len = 100
    ctrate = 10000
    N = len ÷ dt
    dt = len / N
    times = sort(rand(Uniform(0, len), len * ctrate))
    gti = [[0 len];;]
    bins = LinRange(0, len, N + 1)
    counts = fit(Histogram, times, bins).weights
    errs = fill!(similar(counts),1) * sqrt(ctrate)
    bin_times = (@view(bins[1:end-1]) + @view(bins[2:end])) / 2
    segment_size = 20.0
    times2 = sort(rand(Uniform(0, len), len * ctrate))
    counts2 = fit(Histogram, times2, bins).weights
    errs2 = fill!(similar(counts2),1) * sqrt(ctrate)

    @test get_average_ctrate(times, gti, segment_size) == ctrate
    @test get_average_ctrate(bin_times, gti, segment_size; counts=counts) == ctrate

    @testset "test_fts_from_segments_cts_and_events_are_equal" begin
        N = round(Int, segment_size / dt)
        fts_evts = collect(get_flux_iterable_from_segments(times, gti, segment_size, n_bin=N))
        fts_cts = collect(get_flux_iterable_from_segments(
                bin_times, gti, segment_size, fluxes=counts))
        @test fts_evts == fts_cts
    end

    @testset "test_error_on_averaged_cross_spectrum_low_nave" for common_ref in [true, false]
        @test_logs (:warn,r"n_ave is below 30."
        ) error_on_averaged_cross_spectrum([4 + 1.0im], [2], [4], 29, 2, 2; common_ref=common_ref)
    end

    @testset "test_avg_pds_bad_input" begin
        _times = rand(Uniform(0,1000),1)
        out_ev = avg_pds_from_events(_times, gti, segment_size, dt,silent = true)
        @test isnothing(out_ev)
    end

    @testset "test_avg_cs_bad_input" for return_auxil in [true, false]
        _times1 = rand(Uniform(0,1000),1)
        _times2 = rand(Uniform(0,1000),1)
        out_ev = avg_cs_from_events(_times1, _times2, gti,
                                    segment_size, dt, silent = true, return_auxil=return_auxil)
        @test isnothing(out_ev) 
    end

    @testset "test_avg_pds_use_common_mean_similar_stats" for norm in ["frac", "abs", "none", "leahy"]
        out_comm = avg_pds_from_events(
            times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=true,
            silent=true,
            fluxes=nothing,
        ).power
        out = avg_pds_from_events(
            times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=false,
            silent=true,
            fluxes=nothing,
        ).power
        @test std(out_comm)≈std(out) rtol=0.1
    end

    @testset "test_avg_cs_use_common_mean_similar_stats" for 
        norm in ["frac", "abs", "none", "leahy"], 
        return_auxil in [true, false], fullspec in [true,false]
        out_comm = avg_cs_from_events(
            times,
            times2,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=true,
            silent=true,
            fullspec=fullspec,
            return_auxil=return_auxil
        ).power
        out = avg_cs_from_events(
            times,
            times2,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=false,
            silent=true,
            fullspec=fullspec,
            return_auxil=return_auxil
        ).power
        @test std(out_comm)≈std(out) rtol=0.1
    end

    @testset "test_avg_pds_cts_and_events_are_equal" for use_common_mean in [true,false], norm in ["frac", "abs", "none", "leahy"]
        out_ev = avg_pds_from_events(
            times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes=nothing,
        )
        out_ct = avg_pds_from_events(
            bin_times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes=counts,
        )
        compare_tables(out_ev, out_ct)
    end

    @testset "test_avg_pds_cts_and_err_and_events_are_equal" for use_common_mean in [true,false], norm in ["frac", "abs", "none", "leahy"]       
        out_ev = avg_pds_from_events(
            times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes=nothing,
        )
        out_ct = avg_pds_from_events(
            bin_times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes=counts,
            errors=errs,
        )
        # The variance is not _supposed_ to be equal, when we specify errors
        if use_common_mean
            compare_tables(out_ev, out_ct, rtol=0.01, discard=["variance"])
        else
            compare_tables(out_ev, out_ct, rtol=0.1, discard=["variance"])
        end
    end

    @testset "test_avg_cs_cts_and_events_are_equal" for use_common_mean in [true,false], norm in ["frac", "abs", "none", "leahy"]
        out_ev = avg_cs_from_events(
            times,
            times2,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
        )
        out_ct = avg_cs_from_events(
            bin_times,
            bin_times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes1=counts,
            fluxes2=counts2,
        )
        if use_common_mean
            compare_tables(out_ev, out_ct, rtol=0.01)
        else
            compare_tables(out_ev, out_ct, rtol=0.1)
        end
    end

    @testset "test_avg_cs_cts_and_err_and_events_are_equal" for use_common_mean in [true,false], norm in ["frac", "abs", "none", "leahy"]
        out_ev = avg_cs_from_events(
            times,
            times2,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
        )
        out_ct = avg_cs_from_events(
            bin_times,
            bin_times,
            gti,
            segment_size,
            dt,
            norm=norm,
            use_common_mean=use_common_mean,
            silent=true,
            fluxes1=counts,
            fluxes2=counts2,
            errors1=errs,
            errors2=errs2,
        )
        discard = [m for m in propertynames(out_ev) if m == :variance]

        if use_common_mean
            compare_tables(out_ev, out_ct, rtol=0.01, discard=discard)
        else
            compare_tables(out_ev, out_ct, rtol=0.1, discard=discard)
        end
    end
end

@testset "test_norm" begin
    mean = var = 100000
    N = 1000000
    dt = 0.2
    meanrate = mean / dt
    lc = rand(Poisson(mean),N)
    pds = abs2.(fft(lc))
    freq = fftfreq(N, dt)
    good = 2:(N÷2)

    pdsabs = normalize_abs(pds, dt, size(lc,1))
    pdsfrac = normalize_frac(pds, dt, size(lc,1), mean)
    pois_abs = poisson_level("abs", meanrate=meanrate)
    pois_frac = poisson_level("frac", meanrate=meanrate)

    @test Statistics.mean(keepat!(pdsabs,good))≈pois_abs rtol=0.02
    @test Statistics.mean(keepat!(pdsfrac,good))≈pois_frac rtol=0.02

    mean = var = 100000.0
    N = 800000
    dt = 0.2
    df = 1 / (N * dt)
    freq = fftfreq(N, dt)
    good = freq .> 0
    meanrate = mean / dt
    lc = rand(Poisson(mean),N)
    nph = sum(lc)
    pds = keepat!(abs2.(fft(lc)),good)
    lc_bksub = lc .- mean
    pds_bksub = keepat!(abs2.(fft(lc_bksub)),good)
    lc_renorm = lc / mean
    pds_renorm = keepat!(abs2.(fft(lc_renorm)),good)
    lc_renorm_bksub = lc_renorm .- 1
    pds_renorm_bksub = keepat!(abs2.(fft(lc_renorm_bksub)),good)

    @testset "test_leahy_bksub_var_vs_standard" begin
        leahyvar = normalize_leahy_from_variance(pds_bksub, Statistics.var(lc_bksub), N)
        leahy =  2 .* pds ./ sum(lc)
        ratio = Statistics.mean(leahyvar./leahy)
        @test ratio≈1 rtol=0.01
    end
    @testset "test_abs_bksub" begin
        ratio = normalize_abs(pds_bksub, dt, N) ./ normalize_abs(pds, dt, N)
        @test Statistics.mean(ratio)≈1 rtol=0.01
    end
    @testset "test_frac_renorm_constant" begin
        ratio = normalize_frac(pds_renorm, dt, N, 1) ./ normalize_frac(pds, dt, N, mean)
        @test Statistics.mean(ratio)≈1 rtol=0.01
    end
    @testset "test_frac_to_abs_ctratesq" begin
        ratio = (
            normalize_frac(pds, dt, N, mean)
            ./ normalize_abs(pds, dt, N)
            .* meanrate .^ 2
        )
        @test Statistics.mean(ratio)≈1 rtol=0.01
    end
    @testset "test_total_variance" begin
        vdk_total_variance = sum((lc .- mean) .^ 2)
        ratio = Statistics.mean(pds) / vdk_total_variance
        @test Statistics.mean(ratio)≈1 rtol=0.01
    end
    @testset "test_poisson_level_$(norm)" for norm in ["abs", "frac", "leahy","none"]
        pdsnorm = normalize_periodograms(pds, dt, N; mean_flux=mean, n_ph=nph, norm=norm)
        @test Statistics.mean(pdsnorm)≈poisson_level(norm; meanrate=meanrate, n_ph=nph) rtol=0.01
    end

    @testset "test_poisson_level_real_$(norm)" for norm in ["abs", "frac", "leahy","none"]
        pdsnorm = normalize_periodograms(pds, dt, N; mean_flux=mean, n_ph=nph, norm=norm, power_type = "real")
        @test Statistics.mean(pdsnorm)≈poisson_level(norm; meanrate=meanrate, n_ph=nph) rtol=0.01
    end

    @testset "test_poisson_level_absolute_$(norm)" for norm in ["abs", "frac", "leahy","none"]
        pdsnorm = normalize_periodograms(pds, dt, N; mean_flux=mean, n_ph=nph, norm=norm, power_type = "abs")
        @test Statistics.mean(pdsnorm)≈poisson_level(norm; meanrate=meanrate, n_ph=nph) rtol=0.01
    end

    @testset "test_normalize_with_variance" begin
        pdsnorm = normalize_periodograms(pds, dt, N; mean_flux=mean, variance = var, norm="leahy")
        @test Statistics.mean(pdsnorm)≈2 rtol=0.01
    end

    @testset "test_normalize_none" begin
        pdsnorm = normalize_periodograms(pds, dt, N; mean_flux=mean, n_ph=nph, norm="none")
        @test Statistics.mean(pdsnorm)≈Statistics.mean(pds) rtol=0.01
    end
end
