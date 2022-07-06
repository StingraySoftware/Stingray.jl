@testset "test_time_intervals_from_gtis" begin
        start_times, stop_times = time_intervals_from_gtis([[0 400]; [1022 1200];
                                  [1210 1220];], 128)
        @test start_times == [0, 128, 256, 1022]
        @test stop_times == start_times .+ 128
end

@testset "test_time_intervals_from_gtis_frac" begin
        start_times, stop_times = time_intervals_from_gtis([[0 400]; [1022 1200];
                                  [1210 1220];], 128, fraction_step=0.5)
        @test start_times == [0, 64, 128, 192, 256, 1022]
        @test stop_times == start_times .+ 128
end

@testset "test_bin_intervals_from_gtis" begin
    times = range(0.5, 12.5)
    start_bins, stop_bins = bin_intervals_from_gtis([[0 5]; [6 8];], 2, times)

    @test start_bins == [0, 2, 6]
    @test stop_bins == [2, 4, 8]
end

@testset "test_bin_intervals_from_gtis_2" begin
    dt = 0.1
    tstart = 0
    tstop = 100
    times = range(tstart, tstop, step=dt)
    gti = [[tstart - dt/2 tstop - dt/2];]
    start_bins, stop_bins = bin_intervals_from_gtis(gti, 20, times)
    @test start_bins == [0, 200, 400, 600, 800]
end

@testset "test_bin_intervals_from_gtis_frac" begin
    times = range(0.5, 12.5)
    start_bins, stop_bins = bin_intervals_from_gtis([[0 5]; [6 8];], 2, times, fraction_step=0.5)

    @test start_bins == [0, 1, 2, 3, 6]
    @test stop_bins == [2, 3, 4, 5, 8]
end
