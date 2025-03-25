@testset "test_load_gtis" begin
    fname = joinpath(@__DIR__ ,"data","monol_testA.evt")
    @test load_gtis(fname) == [8.0e7 8.0001025e7]
end

@testset "get_gti_length" begin
    @test get_total_gti_length([[0 5]; [6 7];]) == 6
end

@testset "test_check_gtis" begin
    @testset "test_check_gtis_shape" begin
        @test_throws ArgumentError Stingray.check_gtis([[0 1 4]; [0 3 4];]) 
    end
    
    @testset "test_check_gtis_values" begin
        @test_throws ArgumentError Stingray.check_gtis([[0 2]; [1 3];])
    
        @test_throws ArgumentError Stingray.check_gtis([[1 0];])
    end
end

@testset "test_gti_mask" begin

    @testset "test_gti_mask1" begin
        arr = [0, 1, 2, 3, 4, 5, 6]
        gti = [[0 2.1]; [3.9 5];]
        mask, new_gtis = create_gti_mask(arr, gti)
        # NOTE: the time bin has to be fully inside the GTI. That is why the
        # bin at times 0, 2, 4 and 5 are not in.
        @test mask == [0, 1, 0, 0, 0, 0, 0]
    end

    @testset "test_gti_mask_minlen" begin
        arr = [0, 1, 2, 3, 4, 5, 6]
        gti = [[0 2.1]; [3.9 5];]
        mask, new_gtis = create_gti_mask(arr, gti; min_length=2)
        # NOTE: the time bin has to be fully inside the GTI. That is why the
        # bin at times 0, 2, 4 and 5 are not in.
        @test mask == [0, 1, 0, 0, 0, 0, 0]
        @test new_gtis == [0 2.1]
    end
    
    @testset "test_gti_mask_none_longer_than_minlen" begin
        arr = [0, 1, 2, 3, 4, 5, 6]
        gti = [[0 2.1]; [3.9 5];]
        mask = Bool[]
        @test_logs (:warn,r"No GTIs longer than"
            ) mask, _ = create_gti_mask(arr, gti; min_length=10)
        @test all(iszero, mask)
    end
    
    @testset "test_gti_mask_fails_empty_time" begin
        arr = Float64[]
        gti = [[0 2.1]; [3.9 5];]
        @test_throws ArgumentError create_gti_mask(arr, gti)
    end
end

@testset "test_gti_from_condition" begin
    @testset "test_gti_from_condition1" begin
        t = [0, 1, 2, 3, 4, 5, 6]
        condition = [true, true, false, false, true, false, false]
        gti = create_gti_from_condition(t, condition)
        @test gti == [[-0.5 1.5]; [3.5 4.5];]
    end
            
    @testset "test_gti_from_condition2" begin
        t = [0, 1, 2, 3, 4, 5, 6]
        condition = [true, true, true, true, false, true, false]
        gti = create_gti_from_condition(t, condition, safe_interval=[1, 1])
        @test gti == [0.5 2.5]
    end
    
    @testset "test_gti_from_condition_fail" begin
        t = [0, 1, 2, 3]
        condition = [true, true, true]
        @test_throws ArgumentError create_gti_from_condition(t, condition, safe_interval=[1, 1])
    end
end

@testset "test_operations_on_gti" begin
    @testset "test_intersectgti1" begin
        gti1 = [1 4]
        gti2 = [2 5]
        newgti = operations_on_gtis([gti1, gti2], intersect)
        @test newgti == [2 4]
    end
    
    @testset "test_intersectgti2" begin
        gti1 = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
        gti2 = [[2 5]; [6 9]; [11.4 14];]
        newgti = operations_on_gtis([gti1, gti2], intersect)
        @test newgti == [[4.0 5.0]; [7.0 9.0]; [12.2 13.2];]
    end
    
    @testset "test_intersectgti3" begin
        gti1 = [[1 2]; [4 5]; [7 10];]
        newgti = operations_on_gtis([gti1], intersect)
        @test newgti == gti1
    end

    @testset "test_union_gtis_nonoverlapping" begin
        gti0 = [[0 1]; [2 3];]
        gti1 = [[10 11]; [12 13];]
        @test operations_on_gtis([gti0, gti1], union) == [[0 1]; [2 3]; [10 11]; [12 13];]
    end
    
    @testset "test_union_gtis_overlapping" begin
        gti0 = [[0 1]; [2 3]; [4 8];]
        gti1 = [[7 8]; [10 11]; [12 13];]
        @test operations_on_gtis([gti0, gti1], union) == [[0 1]; [2 3]; [4 8]; [10 11]; [12 13];]
    end
end

@testset "test_bti" begin
    @testset "test_bti" begin
        gti = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
        bti = get_btis(gti)
        @test bti == [[2 4]; [5 7]; [10 11]; [11.2 12.2];]
    end    
    
    @testset "test_bti_start_and_stop" begin
        gti = [[1 2]; [4 5]; [7 10]; [11 11.2]; [12.2 13.2];]
        bti = get_btis(gti, 0, 14)
        @test bti == [[0 1]; [2 4]; [5 7]; [10 11]; [11.2 12.2]; [13.2 14];]
    end        
    
    @testset "test_bti_empty_valid" begin
        gti = reshape(Float64[],0,2)
        bti = get_btis(gti, 0, 1)
        @test bti == [0 1]
    end       
    
    @testset "test_bti_fail" begin
        gti = reshape(Float64[],0,2)
        @test_throws ArgumentError get_btis(gti)
    end  
end

@testset "test_time_intervals_from_gtis" begin
    @testset "test_time_intervals_from_gtis1" begin
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
end

@testset "test_bin_intervals_from_gtis" begin
    @testset "test_bin_intervals_from_gtis1" begin
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
end

@testset "GTI Interface" begin
    times = collect(0.0:0.1:1.0) 
    energies = rand(length(times))
    el = EventList("test.evt", times, energies, DictMetadata([Dict()]))

    lc = create_lightcurve(el, 0.1)
    @test length(lc.timebins) == 10
    @test length(lc.counts) == 10
    @test length(lc.count_error) == 10
    
    gtis = [0.05 0.25; 0.35 0.55; 0.65 0.85]  

    @testset "apply_gtis LightCurve" begin
        filtered = apply_gtis(lc, gtis)
        @test length(filtered) == size(gtis, 1)
        
        for (i, f) in enumerate(filtered)
            @test length(f.timebins) == length(f.counts) == length(f.count_error)
            @test all(f.timebins .>= gtis[i,1])
            @test all(f.timebins .<= gtis[i,2])
        end
    end
    
    @testset "fill_bad_time_intervals! LightCurve" begin
        test_lc = deepcopy(lc)
        fill_bad_time_intervals!(test_lc, gtis, fill_value=-1.0)
        
        @test any(test_lc.counts .== -1.0)
        @test all(test_lc.count_error[test_lc.counts .== -1.0] .== 0.0)
        
        for i in 1:size(gtis,1)
            in_gti = (test_lc.timebins .>= gtis[i,1]) .& (test_lc.timebins .<= gtis[i,2])
            @test all(test_lc.counts[in_gti] .!= -1.0)
        end
    end
end
