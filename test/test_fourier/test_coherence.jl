function test_data()
    fid = h5open(joinpath(@__DIR__, "..", "data", "sample_variable_lc.h5"), "r")
    data = (fid["dataset1"][1:10000]) .* 1000.0
    close(fid)
    dt = 0.01
    N = length(data)
    mean = Statistics.mean(data)
    meanrate = mean / dt
    freq = fftfreq(length(data), 1/dt)
    good = 0 .< freq .< 0.1
    return (data=data, dt=dt, N=N, mean=mean, meanrate=meanrate, good=good)
end

# Helper function to generate processed data for tests
function generate_test_data(data; same_poisson=true)
    data1 = rand.(Poisson.(data.data))
    data2 = same_poisson ? data1 : rand.(Poisson.(data.data))
    
    ft1 = fft(data1)
    ft2 = fft(data2)
    keepat!(ft1, data.good)
    keepat!(ft2, data.good)
    cross = normalize_periodograms(
        ft1 .* conj(ft2), data.dt, data.N, 
        mean_flux = data.mean, norm="abs", power_type="all")
    pds1 = normalize_periodograms(
        ft1 .* conj(ft1), data.dt, data.N, 
        mean_flux = data.mean, norm="abs", power_type="real")
    pds2 = normalize_periodograms(
        ft2 .* conj(ft2), data.dt, data.N, 
        mean_flux = data.mean, norm="abs", power_type="real")
    p1noise = poisson_level("abs", meanrate=data.meanrate)
    p2noise = poisson_level("abs", meanrate=data.meanrate)
    
    return (cross=cross, pds1=pds1, pds2=pds2, p1noise=p1noise, p2noise=p2noise)
end
data = test_data()

# test_coherence
let
    test_data = generate_test_data(data, same_poisson=false)
end

# test_intrinsic_coherence
let
    test_data = generate_test_data(data, same_poisson=true)
    
    coh = estimate_intrinsic_coherence.(
        test_data.cross, test_data.pds1, test_data.pds2, 
        test_data.p1noise, test_data.p2noise, data.N)
    @test all(x -> isapprox(x, 1; atol=0.001), coh)
end

# test_raw_high_coherence  
let
    test_data = generate_test_data(data, same_poisson=true)
    
    coh = raw_coherence.(test_data.cross, test_data.pds1, test_data.pds2, 
                        test_data.p1noise, test_data.p2noise, data.N)
    @test all(x -> isapprox(x, 1; atol=0.001), coh)
end

# test_raw_low_coherence
let
    test_data = generate_test_data(data, same_poisson=true)
    
    nbins = 2
    C, P1, P2 = @view(test_data.cross[1:nbins]), @view(test_data.pds1[1:nbins]), @view(test_data.pds2[1:nbins])
    bsq = bias_term.(P1, P2, test_data.p1noise, test_data.p2noise, data.N)
    # must be lower than bsq!
    low_coh_cross = @. rand(Normal(bsq^0.5 / 10, bsq^0.5 / 100)) + 0.0im
    coh = raw_coherence.(low_coh_cross, P1, P2, test_data.p1noise, test_data.p2noise, data.N)
    @test all(iszero, coh)
end

# test_raw_high_bias (this one was already different, so keeping it separate)
let
    C = [12986.0 + 8694.0im]
    P1 = [476156.0]
    P2 = [482751.0]
    P1noise = 495955
    P2noise = 494967

    coh = raw_coherence.(C, P1, P2, P1noise, P2noise, 499; intrinsic_coherence=1)
    coh_sngl = raw_coherence(C[1], P1[1], P2[1], P1noise, P2noise, 499; intrinsic_coherence=1)
    @test coh == real(C .* conj(C)) ./ (P1 .* P2)
    @test coh_sngl == real(C .* conj(C))[1] / (P1[1] * P2[1])
end