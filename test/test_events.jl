@testset "test_events" begin
    time = [0.5, 1.5, 2.5, 3.5]
    counts = [3000, 2000, 2200, 3600]
    counts_flat = [3000, 3000, 3000, 3000]
    spectrum = [[1, 2, 3, 4, 5, 6], [1000, 2040, 1000, 3000, 4020, 2070]]
    gti = [0 4]

    @testset "test_initialize_eventList" begin
        times = sort(rand(Uniform(1e8,1e8+100),101))
        ev = EventList(time=times, mjdref=54600)
        @test ev.time ≈ times atol = 1e-15
        @test ev.mjdref == 54600 
    end

    @testset "test_to_and_from_lightcurve" begin
        @testset "test_from_lc" begin
            lc = LightCurve(time=[0.5, 1.5, 2.5], counts=[2, 1, 2])
            ev = from_lc(lc)
            @test ev.time == [0.5, 0.5, 1.5, 2.5, 2.5]
        end

        @testset "test_to_lc" begin
            ev = EventList(time=time, gti=gti)
            lc = to_lc(ev, 1)
            @test lc.time == [0.5, 1.5, 2.5, 3.5]
            @test lc.gti == gti
        end
    end

    @testset "test_join_events" begin
        @testset "test_join_different_dt" begin
            ev = EventList(time=[10, 20, 30], dt=1)
            ev_other = EventList(time=[40, 50, 60], dt=3)
            ev_new = Stingray.join(ev,ev_other)
            @test ev_new.dt == 3
        end
        @testset "test_join_different_instr" begin
            ev = EventList(time=[10, 20, 30], instr="fpma")
            ev_other = EventList(time=[40, 50, 60], instr="fpmb")
            ev_new = Stingray.join(ev,ev_other)
            @test ev_new.instr == "fpma,fpmb"
        end
        @testset "test_join_without_energy" begin
            ev = EventList(time=[1, 2, 3], energy=[3, 3, 3])
            ev_other = EventList(time=[4, 5])
            ev_new = Stingray.join(ev, ev_other)
            @test ev_new.energy == [3, 3, 3, 0, 0]
        end
        @testset "test_join_without_pi" begin
            ev = EventList(time=[1, 2, 3], PI=[3, 3, 3])
            ev_other = EventList(time=[4, 5])
            ev_new = Stingray.join(ev, ev_other)
            @test ev_new.PI == [3, 3, 3, 0, 0]
        end
        @testset "test_join_with_gti_none" begin
            ev = EventList(time=[1, 2, 3])
            ev_other = EventList(time=[4, 5], gti=[3.5 5.5])
            ev_new = Stingray.join(ev, ev_other)
            @test ev_new.gti == [[1 3]; [3.5 5.5];]

            ev = EventList(time=[1, 2, 3], gti=[0.5 3.5])
            ev_other = EventList(time=[4, 5])
            ev_new = Stingray.join(ev,ev_other)
            @test ev_new.gti == [[0.5 3.5]; [4 5];]

            ev = EventList(time=[1, 2, 3])
            ev_other = EventList(time=[4, 5])
            ev_new = Stingray.join(ev, ev_other)
            @test isempty(ev_new.gti)
        end
        @testset "test_non_overlapping_join" begin
            ev = EventList(time=[1, 1, 2, 3, 4],
                        energy=[3, 4, 7, 4, 3], gti=[[1 2]; [3 4];])
            ev_other = EventList(time=[5, 6, 6, 7, 10],
                                energy=[4, 3, 8, 1, 2], gti=[6 7])
            ev_new = Stingray.join(ev, ev_other)

            @test ev_new.time == [1, 1, 2, 3, 4, 5, 6, 6, 7, 10]
            @test ev_new.energy == [3, 4, 7, 4, 3, 4, 3, 8, 1, 2]
            @test ev_new.gti == [[1 2]; [3 4]; [6 7];]
        end
        @testset "test_overlapping_join" begin
            ev = EventList(time=[1, 1, 10, 6, 5],
                        energy=[10, 6, 3, 11, 2], gti=[[1 3]; [5 6];])
            ev_other = EventList(time=[5, 7, 6, 6, 10],
                                energy=[2, 3, 8, 1, 2], gti=[[5 7]; [8 10];])
            ev_new = Stingray.join(ev, ev_other)

            @test ev_new.time == [1, 1, 5, 5, 6, 6, 6, 7, 10, 10]
            @test ev_new.energy == [10, 6, 2, 2, 11, 8, 1, 3, 3, 2]
            @test ev_new.gti == [5 6]
        end
    end
    @testset "test_sort" begin
        ev = EventList(time=[2, 4, 1, 4, 3], energy=[3, 4, 7, 4, 3])
        @testset "test_sort_not_inplace" begin
            new_ev = Stingray.sort(ev);
            @test new_ev.time == [1, 2, 3, 4, 4]
            @test new_ev.energy == [7, 3, 3, 4, 4]
            @test ev.time == [2, 4, 1, 4, 3]
        end
        @testset "test_sort_inplace" begin
            Stingray.sort!(ev);
            @test ev.time == [1, 2, 3, 4, 4]
            @test ev.energy == [7, 3, 3, 4, 4] 
        end
    end
    @testset "Input/Output" begin
        @testset "test_fits_with_standard_file" begin
            fname = joinpath(@__DIR__ ,"data","monol_testA.evt")
            ev = read(EventList, fname, "fits")
            @test ev.mjdref == 55197.00076601852
        end

        @testset "test_io_with_fits" begin
            ev = EventList(time=time, mjdref=54000)
            write(ev, "ev.fits", "fits")
            new_ev = read(EventList, "ev.fits", "fits")
            @test new_ev.time == time
            rm("ev.fits")
        end
    end
end

