using Test
using Stingray

@testset "Powerspectrum" begin

    @testset "Direct struct initialization" begin
        freq = [1.0, 2.0, 3.0]
        power = [1.5, 2.0, 1.8]
        power_err = [0.15, 0.2, 0.18]
        unnorm = [150.0, 200.0, 180.0]
        unnorm_err = [15.0, 20.0, 18.0]
        test_gti = [0.0 10.0]

        ps = Powerspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            1.0, 0.0001, 100000, 1, 1,  # df, dt, n, m, k
            1000.0,  # nphots
            "leahy",  # norm
            test_gti, 10.0,  # gti, segment_size
            nothing,  # variance
            "poisson", "powerspectrum"  # err_dist, type
        )

        @test ps isa Powerspectrum{Float64}
        @test ps isa AbstractPowerspectrum
        @test ps isa AbstractCrossspectrum  # subtype check
        @test length(ps.freq) == 3
        @test ps.norm == "leahy"
        @test ps.type == "powerspectrum"
        @test ps.m == 1
        @test ps.k == 1
        @test isnothing(ps.variance)
        @test ps.err_dist == "poisson"
        @test eltype(ps.power) <: Real
    end

end
