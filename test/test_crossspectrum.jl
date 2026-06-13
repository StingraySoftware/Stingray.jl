using Test
using Stingray

@testset "Crossspectrum" begin

    @testset "Direct struct initialization" begin
        freq = [1.0, 2.0, 3.0]
        power = Complex{Float64}.([1.0+0.5im, 2.0+1.0im, 3.0+0.2im])
        power_err = Complex{Float64}.([0.1+0.05im, 0.2+0.1im, 0.3+0.02im])
        unnorm = Complex{Float64}.([10.0+5.0im, 20.0+10.0im, 30.0+2.0im])
        unnorm_err = Complex{Float64}.([1.0+0.5im, 2.0+1.0im, 3.0+0.2im])
        test_gti = [0.0 10.0]

        cs = Crossspectrum{Float64}(
            freq, power, power_err, unnorm, unnorm_err,
            nothing, nothing,  # pds1, pds2
            1.0, 0.001, 10000, 1, 1,  # df, dt, n, m, k
            100.0, 120.0, 80.0,  # nphots, nphots1, nphots2
            12.0, 8.0,  # countrate1, countrate2
            "frac", "all", false,  # norm, power_type, fullspec
            test_gti, 10.0,  # gti, segment_size
            "poisson", "crossspectrum"  # err_dist, type
        )

        @test cs isa Crossspectrum{Float64}
        @test cs isa AbstractCrossspectrum
        @test length(cs.freq) == 3
        @test cs.norm == "frac"
        @test cs.type == "crossspectrum"
        @test cs.m == 1
        @test cs.k == 1
        @test cs.fullspec == false
        @test cs.err_dist == "poisson"
        @test eltype(cs.power) <: Complex
    end

end
