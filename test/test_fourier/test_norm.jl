const TEST_MEAN = TEST_VAR = 100000.0
const TEST_N = 800000
const TEST_DT = 0.2

function setup_test_data(; N=TEST_N, dt=TEST_DT, mean=TEST_MEAN)
    freq = fftfreq(N, dt)
    good = freq .> 0
    meanrate = mean / dt
    lc = rand(Poisson(mean), N)
    nph = sum(lc)
    pds = keepat!(abs2.(fft(lc)), good)
    return (
        freq=freq, good=good, meanrate=meanrate, 
        lc=lc, nph=nph, pds=pds, N=N, dt=dt, mean=mean
    )
end

# Helper function for Poisson level tests
function test_poisson_level(norm::String, power_type::String="all")
    data = setup_test_data()
    pdsnorm = normalize_periodograms(
        data.pds, data.dt, data.N; 
        mean_flux=data.mean, n_ph=data.nph, 
        norm=norm, power_type=power_type
    )
    expected = poisson_level(norm; meanrate=data.meanrate, n_ph=data.nph)
    @test Statistics.mean(pdsnorm) ≈ expected rtol=0.01
end

# test_positive_fft_bins
let
    freq = fftfreq(11)
    goodbins = positive_fft_bins(11)
    @test filter(x -> x > 0, freq) == freq[goodbins]
end

# test_norm_setup_and_basic_tests
let
    data = setup_test_data(N=1000000)
    pdsabs = normalize_abs(data.pds, data.dt, size(data.lc, 1))
    pdsfrac = normalize_frac(data.pds, data.dt, size(data.lc, 1), data.mean)
    pois_abs = poisson_level("abs", meanrate=data.meanrate)
    pois_frac = poisson_level("frac", meanrate=data.meanrate)

    @test Statistics.mean(pdsabs) ≈ pois_abs rtol=0.02
    @test Statistics.mean(pdsfrac) ≈ pois_frac rtol=0.02
end

# test_leahy_bksub_var_vs_standard
let
    data = setup_test_data()
    lc_bksub = data.lc .- data.mean
    pds_bksub = keepat!(abs2.(fft(lc_bksub)), data.good)
    leahyvar = normalize_leahy_from_variance(pds_bksub, Statistics.var(lc_bksub), data.N)
    leahy = 2 .* data.pds ./ sum(data.lc)
    ratio = Statistics.mean(leahyvar ./ leahy)
    @test ratio ≈ 1 rtol=0.01
end

# test_abs_bksub
let
    data = setup_test_data()
    lc_bksub = data.lc .- data.mean
    pds_bksub = keepat!(abs2.(fft(lc_bksub)), data.good)
    ratio = normalize_abs(pds_bksub, data.dt, data.N) ./ normalize_abs(data.pds, data.dt, data.N)
    @test Statistics.mean(ratio) ≈ 1 rtol=0.01
end
# test_frac_renorm_constant
let
    data = setup_test_data()
    lc_renorm = data.lc / data.mean
    pds_renorm = keepat!(abs2.(fft(lc_renorm)), data.good)
    ratio = normalize_frac(pds_renorm, data.dt, data.N, 1) ./ normalize_frac(data.pds, data.dt, data.N, data.mean)
    @test Statistics.mean(ratio) ≈ 1 rtol=0.01
end

# test_frac_to_abs_ctratesq
let
    data = setup_test_data()
    
    ratio = (
        normalize_frac(data.pds, data.dt, data.N, data.mean)
        ./ normalize_abs(data.pds, data.dt, data.N)
        .* data.meanrate .^ 2
    )
    @test Statistics.mean(ratio) ≈ 1 rtol=0.01
end

# test_total_variance
let
    data = setup_test_data()
    
    vdk_total_variance = sum((data.lc .- data.mean) .^ 2)
    ratio = Statistics.mean(data.pds) / vdk_total_variance
    @test Statistics.mean(ratio) ≈ 1 rtol=0.01
end

# All Poisson level tests - much more concise now!
@testset "Poisson Level Tests" begin
    test_poisson_level("abs")
    test_poisson_level("frac") 
    test_poisson_level("leahy")
    test_poisson_level("none")
    
    # Real power type tests
    test_poisson_level("abs", "real")
    test_poisson_level("frac", "real")
    test_poisson_level("leahy", "real") 
    test_poisson_level("none", "real")
    
    # Absolute power type tests
    test_poisson_level("abs", "abs")
    test_poisson_level("frac", "abs")
    test_poisson_level("leahy", "abs")
    test_poisson_level("none", "abs")
end

# test_normalize_with_variance
let
    data = setup_test_data()
    pdsnorm = normalize_periodograms(data.pds, data.dt, data.N; mean_flux=data.mean, variance=TEST_VAR, norm="leahy")
    @test Statistics.mean(pdsnorm) ≈ 2 rtol=0.01
end

# test_normalize_none
let
    data = setup_test_data()
    pdsnorm = normalize_periodograms(data.pds, data.dt, data.N; mean_flux=data.mean, n_ph=data.nph, norm="none")
    @test Statistics.mean(pdsnorm) ≈ Statistics.mean(data.pds) rtol=0.01
end