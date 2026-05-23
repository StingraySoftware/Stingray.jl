using Random
using Distributions

rng = Random.Xoshiro(1259723)

@testset "Power Colors" begin

    freq = collect(0.0001:0.00001:17)
    pow_spec = 1.0 ./ freq

    pc0, pc0e, pc1, pc1e = power_color(freq, pow_spec)
    lpc0, _, lpc1, _ = power_color(freq, pow_spec; return_log = true)

    @testset "power_color: 1/f spectrum gives pc ≈ 1" begin
        @test isapprox(pc0, 1.0; atol = 1e-3)
        @test isapprox(pc1, 1.0; atol = 1e-3)
    end

    @testset "power_color: log10(pc) ≈ 0 for 1/f" begin
        @test isapprox(lpc0, 0.0; atol = 0.001)
        @test isapprox(lpc1, 0.0; atol = 0.001)
    end

    @testset "power_color: bad frequency edges" begin
        good = freq .> (1 / 255)
        @test_throws ArgumentError power_color(freq[good], pow_spec[good])

        good = freq .< 15
        @test_throws ArgumentError power_color(freq[good], pow_spec[good])

        @test_throws ArgumentError power_color(freq, pow_spec; freq_edges = [1.0])
        @test_throws ArgumentError power_color(freq, pow_spec;
            freq_edges = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    end

    @testset "power_color: bad excluded interval format" begin
        for fte in [
            [1.0, 1.1, 3.0],
            [[1.0, 1.1, 3.0]],
        ]
            @test_throws ArgumentError power_color(freq, pow_spec;
                freqs_to_exclude = fte)
        end
    end

    @testset "power_color: excluding small interval" begin
        pc0_ex, _, pc1_ex, _ = power_color(freq, pow_spec;
            freqs_to_exclude = [(1.0, 1.1)])
        @test isapprox(pc0_ex, 1.0; atol = 0.001)
        @test isapprox(pc1_ex, 1.0; atol = 0.001)
    end

    @testset "power_color: error propagation" begin
        pc0_h, pc0_err_h, pc1_h, pc1_err_h = power_color(freq, pow_spec;
            power_err = pow_spec ./ 2)
        pc0_f, pc0_err_f, pc1_f, pc1_err_f = power_color(freq, pow_spec;
            power_err = pow_spec)

        @test isapprox(pc0_h, 1.0; atol = 0.001)
        @test isapprox(pc1_h, 1.0; atol = 0.001)
        @test isapprox(pc0_f, 1.0; atol = 0.001)
        @test isapprox(pc1_f, 1.0; atol = 0.001)
        @test isapprox(pc0_err_f / pc0_err_h, 2.0; atol = 0.001)
        @test isapprox(pc1_err_f / pc1_err_h, 2.0; atol = 0.001)
    end

    @testset "hue: round-trip angle recovery" begin
        center = (4.51920, 0.453724)
        log_center = log10.([center[1], center[2]])

        for angle_deg in 0:20:380
            angle = deg2rad(angle_deg)
            factor = rand(rng, Uniform(0.1, 10))
            x = factor * cos(3 / 4 * π - angle) + log_center[1]
            y = factor * sin(3 / 4 * π - angle) + log_center[2]

            hue = hue_from_power_color(10^x, 10^y; center = center)

            c2 = (sin(hue) - sin(angle))^2 + (cos(hue) - cos(angle))^2
            angle_diff = acos(clamp((2.0 - c2) / 2.0, -1.0, 1.0))

            @test isapprox(angle_diff, 0.0; atol = 0.001)
        end
    end

    @testset "integrate_power_in_frequency_range" begin
        n = 1000
        f = collect(range(0.1, 10.0; length = n))
        p = fill(2.0, n)
        frange = [1.0, 5.0]
        val, _ = integrate_power_in_frequency_range(f, p, frange)
        @test isapprox(val, 2.0 * (5.0 - 1.0); rtol = 0.01)

        val_pois, _ = integrate_power_in_frequency_range(f, p, frange;
            poisson_power = 1.0)
        @test isapprox(val_pois, 1.0 * (5.0 - 1.0); rtol = 0.01)
    end

end
