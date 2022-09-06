@testset "test_lightcurve" begin
    @testset "rebin_lightcurve" begin
        dt = 0.0001220703125
        n = 1384132
        mean_counts = 2.0
        times = range(dt/2, dt/2 + n*dt, step = dt)
        counts = zero(times) .+ mean_counts
        lc = LightCurve(time=times, counts=counts)

        @testset "test_rebin_even" begin
            dt_new = 2.0
            lc_binned = rebin(lc, dt_new)
            @test lc_binned.dt == dt_new
            counts_test = zero(lc_binned.time) .+ lc.counts[1]*dt_new/lc.dt
            @test lc_binned.counts == counts_test
        end

        @testset "test_rebin_odd" begin
            dt_new = 1.5
            lc_binned = rebin(lc, dt_new)
            @test lc_binned.dt == dt_new
            counts_test = zero(lc_binned.time) .+ lc.counts[1]*dt_new/lc.dt
            @test lc_binned.counts == counts_test
        end

        @testset "test_rebin_several" for dt in [2,3,pi,5]
            lc_binned = rebin(lc, dt)
            @test length(lc_binned.time) == length(lc_binned.counts)
        end

        @testset "test_rebin_with_gtis" begin
            times = collect(range(0, 99.9, step=0.1))

            counts = round.(rand(Normal(100, 0.1), length(times)))
            gti = [[0 40]; [60 100];]

            good, gti = create_gti_mask(times, gti)

            keepat!(times, good)
            keepat!(counts, good)

            lc = LightCurve(time=times,counts=counts, gti=gti, skip_checks=true, dt=0.1)

            lc_rebin = rebin(lc, 1.0)

            @test (lc_rebin.time[40] - lc_rebin.time[39]) >= 1.0
        end
    end
end
