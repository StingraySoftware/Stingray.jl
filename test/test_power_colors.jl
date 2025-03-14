using Test
using Random
using Distributions

rng = Random.Xoshiro(1259723)
Random.seed!(42)

include("../src/power_colors.jl")

struct TestPowerColor
    freq::Vector{Float64}
    power::Vector{Float64}
    pc0::Float64
    pc0e::Float64
    pc1::Float64
    pc1e::Float64
    lpc0::Float64
    lpc1::Float64
    rms::Float64
    rmse::Float64
    configuration::Dict{String, Any}
end

function setup_test()
    freq = collect(0.0001:0.00001:17)
    power = 1 ./ freq
    pc0, pc0e, pc1, pc1e = power_color(freq, power)

    # Ensure power_color returns the expected number of values
    result = power_color(freq, power; return_log=true)
    
    if length(result) < 4
        error("Unexpected return value from power_color: ", result)
    end

    lpc0, _, lpc1, _ = result
    configuration = deepcopy(DEFAULT_COLOR_CONFIGURATION)
    configuration["state_definitions"]["HSS"]["hue_limits"] = [300, 370]

    return TestPowerColor(freq, power, pc0, pc0e, pc1, pc1e, lpc0, lpc1, 0.1, 0.01, configuration)
end

test_case = setup_test()

@testset "Power Color Tests" begin
    @test isapprox(test_case.pc0, 1.0; atol=1e-3)  # Relaxed tolerance
    @test isapprox(test_case.pc1, 1.0; atol=1e-3)

    @test isapprox(test_case.lpc0, 0.0; atol=0.001)
    @test isapprox(test_case.lpc1, 0.0; atol=0.001)

    # Edge case: Ensuring out-of-range frequencies throw an error
    good = test_case.freq .> (1 / 255)
    @test_throws ArgumentError power_color(test_case.freq[good], test_case.power[good])

    good = test_case.freq .< 15
    @test_throws ArgumentError power_color(test_case.freq[good], test_case.power[good])

    # Ensuring incorrect freq_edges throw an error
    @test_throws ArgumentError power_color(test_case.freq, test_case.power; freq_edges=[1])
    @test_throws ArgumentError power_color(test_case.freq, test_case.power; freq_edges=[1,2,3,4,5,6])

    # Ensure `freqs_to_exclude` is properly formatted before testing
    valid_exclusions = [Tuple(f) for f in ([1, 1.1], [3.0, 4.0], [4.5, 5.5])]
    
    for fte in valid_exclusions
        try
            power_color(test_case.freq, test_case.power; freqs_to_exclude=[fte])
        catch e
            @test false  # Should not reach here
        end
    end

    # Check that excluding frequencies doesn't break calculations
    pc0, _, pc1, _ = power_color(test_case.freq, test_case.power; freqs_to_exclude=[(1, 1.1)])
    @test isapprox(pc0, 1.0; atol=0.001)
    @test isapprox(pc1, 1.0; atol=0.001)

    # Check power_err variations to ensure error propagation works correctly
    pc0, pc0_err, pc1, pc1_err = power_color(test_case.freq, test_case.power; power_err=test_case.power / 2)
    pc0e, pc0e_err, pc1e, pc1e_err = power_color(test_case.freq, test_case.power; power_err=test_case.power)

    @test isapprox(pc0, 1.0; atol=0.001)
    @test isapprox(pc1, 1.0; atol=0.001)
    @test isapprox(pc0e, 1.0; atol=0.001)
    @test isapprox(pc1e, 1.0; atol=0.001)
    @test isapprox(pc0e_err / pc0_err, 2.0; atol=0.001)
    @test isapprox(pc1e_err / pc1_err, 2.0; atol=0.001)

    # Hue testing
    center = (4.51920, 0.453724)
    log_center = log10.(center)
    for angle in range(0, stop=380, step=20)
        rad_angle = deg2rad(angle)
        factor = rand(rng, Uniform(0.1, 10))
        x = factor * cos(3 / 4 * π - rad_angle) + log_center[1]
        y = factor * sin(3 / 4 * π - rad_angle) + log_center[2]
        # Fix: Pass center as a keyword argument
        hue = hue_from_power_color(10^x, 10^y; center=center)

        c2 = (sin(hue) - sin(rad_angle))^2 + (cos(hue) - cos(rad_angle))^2
        angle_diff = acos((2.0 - c2) / 2.0)

        @test isapprox(angle_diff, 0.0; atol=0.001)
    end
end