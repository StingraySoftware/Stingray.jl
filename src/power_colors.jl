"""
    DEFAULT_COLOR_CONFIGURATION

Default power color configuration dictionary, containing:
- `"center"`: The center coordinates `[4.51920, 0.453724]` in the power color diagram.
- `"ref_angle"`: The reference angle `3π/4`.
- `"state_definitions"`: Hue limits and colors for HSS, LHS, HIMS, SIMS states.
- `"rms_spans"`: Fractional rms ranges as a function of hue angle.

See Heil et al. 2015, MNRAS, 448, 3348.
"""
const DEFAULT_COLOR_CONFIGURATION = Dict{String, Any}(
    "center" => [4.51920, 0.453724],
    "ref_angle" => 3 * π / 4,
    "state_definitions" => Dict{String, Any}(
        "HSS"  => Dict{String, Any}("hue_limits" => [300, 360], "color" => "red"),
        "LHS"  => Dict{String, Any}("hue_limits" => [-20, 140], "color" => "blue"),
        "HIMS" => Dict{String, Any}("hue_limits" => [140, 220], "color" => "green"),
        "SIMS" => Dict{String, Any}("hue_limits" => [220, 300], "color" => "yellow"),
    ),
    "rms_spans" => Dict{Int, Vector{Float64}}(
        -20  => [0.3, 0.7],
          0  => [0.3, 0.7],
         10  => [0.3, 0.6],
         40  => [0.25, 0.4],
        100  => [0.25, 0.35],
        150  => [0.2, 0.3],
        170  => [0.0, 0.3],
        200  => [0.0, 0.15],
        370  => [0.0, 0.15],
    ),
)

"""
    power_color(frequency, power; kwargs...)

Calculate two power colors from a power spectrum.

Power colors are defined as the ratio of the integrated power (variance) in
adjacent frequency ranges. Given five frequency edges `[f0, f1, f2, f3, f4]`:

- `PC0 = Var([f0, f1]) / Var([f2, f3])`
- `PC1 = Var([f1, f2]) / Var([f3, f4])`

# Arguments
- `frequency::AbstractVector{<:Real}`: Frequencies of the power spectrum.
- `power::AbstractVector{<:Real}`: Power at each frequency.

# Keyword Arguments
- `power_err::Union{Nothing, AbstractVector}=nothing`: Power error bars.
  If `nothing`, defaults to `power / √m`.
- `freq_edges::Vector{Float64}=[1/256, 1/32, 0.25, 2.0, 16.0]`: Five edges
  defining four frequency intervals.
- `df::Union{Nothing, Real}=nothing`: Frequency resolution. If `nothing`,
  computed from `median(diff(frequency))`.
- `m::Int=1`: Number of averaged segments/bins.
- `freqs_to_exclude::Union{Nothing, Vector}=nothing`: Frequency ranges to exclude
  (e.g., QPO frequencies). Format: `[(f0, f1), (f2, f3), ...]` or `(f0, f1)`.
- `poisson_power::Union{Real, AbstractVector}=0`: Poisson noise level.
- `return_log::Bool=false`: If `true`, return `log10` of power colors.

# Returns
`(PC0, PC0_err, PC1, PC1_err)` — the two power colors and their errors.

# References
Heil et al. 2015, MNRAS, 448, 3348.
"""
function power_color(
    frequency::AbstractVector{<:Real},
    power::AbstractVector{<:Real};
    power_err::Union{Nothing, AbstractVector} = nothing,
    freq_edges::Vector{Float64} = [1 / 256, 1 / 32, 0.25, 2.0, 16.0],
    df::Union{Nothing, Real} = nothing,
    m::Int = 1,
    freqs_to_exclude::Union{Nothing, Vector} = nothing,
    poisson_power::Union{Real, AbstractVector} = 0,
    return_log::Bool = false,
)
    if length(freq_edges) != 5
        throw(ArgumentError("freq_edges must have 5 elements"))
    end

    freq = collect(float.(frequency))
    pow = collect(float.(power))

    df_val = isnothing(df) ? median(diff(freq)) : float(df)

    input_frequency_low_edges = freq .- df_val / 2
    input_frequency_high_edges = freq .+ df_val / 2

    if minimum(freq_edges) < input_frequency_low_edges[1]
        throw(ArgumentError("The minimum frequency is larger than the first frequency edge"))
    end
    if maximum(freq_edges) > input_frequency_high_edges[end]
        throw(ArgumentError("The maximum frequency is lower than the last frequency edge"))
    end

    pow_err = isnothing(power_err) ? pow ./ sqrt(m) : collect(float.(power_err))

    if !isnothing(freqs_to_exclude)
        exclusions = if freqs_to_exclude isa Tuple{<:Real, <:Real}
            [(Float64(freqs_to_exclude[1]), Float64(freqs_to_exclude[2]))]
        elseif length(freqs_to_exclude) == 2 && freqs_to_exclude[1] isa Real
            [(Float64(freqs_to_exclude[1]), Float64(freqs_to_exclude[2]))]
        elseif all(x -> (x isa Tuple || x isa AbstractVector) && length(x) == 2, freqs_to_exclude)
            [(Float64(x[1]), Float64(x[2])) for x in freqs_to_exclude]
        else
            throw(ArgumentError("freqs_to_exclude must be of format [[f0, f1], [f2, f3], ...]"))
        end

        for (f0, f1) in exclusions
            excl_mask = @. (input_frequency_low_edges > f0) & (input_frequency_high_edges < f1)
            idx0 = clamp(searchsortedfirst(freq, f0), 1, length(freq))
            idx1 = clamp(searchsortedfirst(freq, f1), 1, length(freq))
            pow[excl_mask] .= (pow[idx0] + pow[idx1]) / 2
        end
    end

    var00, var00_err = integrate_power_in_frequency_range(
        freq, pow, freq_edges[1:2];
        power_err = pow_err, df = df_val, m = m, poisson_power = poisson_power,
    )
    var01, var01_err = integrate_power_in_frequency_range(
        freq, pow, freq_edges[3:4];
        power_err = pow_err, df = df_val, m = m, poisson_power = poisson_power,
    )
    var10, var10_err = integrate_power_in_frequency_range(
        freq, pow, freq_edges[2:3];
        power_err = pow_err, df = df_val, m = m, poisson_power = poisson_power,
    )
    var11, var11_err = integrate_power_in_frequency_range(
        freq, pow, freq_edges[4:5];
        power_err = pow_err, df = df_val, m = m, poisson_power = poisson_power,
    )

    pc0 = var00 / var01
    pc1 = var10 / var11
    pc0_err = pc0 * (var00_err / var00 + var01_err / var01)
    pc1_err = pc1 * (var10_err / var10 + var11_err / var11)

    if return_log
        pc0_err = (1 / pc0) * pc0_err
        pc1_err = (1 / pc1) * pc1_err
        pc0 = log10(pc0)
        pc1 = log10(pc1)
    end

    return pc0, pc0_err, pc1, pc1_err
end

"""
    hue_from_power_color(pc0, pc1; center=(4.51920, 0.453724))

Measure the angle (hue) of a point in the log-power color diagram with respect
to the center.

The inputs `pc0` and `pc1` are in **linear** (not log) scale — the function takes
the log10 internally.

# Arguments
- `pc0`: The (linear, not log!) first power color. Can be scalar or array.
- `pc1`: The (linear, not log!) second power color. Can be scalar or array.

# Keyword Arguments
- `center`: The center coordinates `(PC0_center, PC1_center)` in linear scale.
  Default: `(4.51920, 0.453724)`.

# Returns
- `hue`: Angle(s) in radians.

# References
Heil et al. 2015, MNRAS, 448, 3348.
"""
function hue_from_power_color(pc0, pc1; center = (4.51920, 0.453724))
    log_pc0 = log10.(pc0)
    log_pc1 = log10.(pc1)
    log_center = log10.(center)
    return hue_from_logpower_color(log_pc0, log_pc1; center = (log_center[1], log_center[2]))
end

"""
    hue_from_logpower_color(log10pc0, log10pc1; center=(log10(4.51920), log10(0.453724)))

Measure the angle (hue) of a point in the log-power color diagram with respect
to the center.

Angles are measured in radians, **in the clockwise direction**, with respect to
a line oriented at −45° with respect to the horizontal axis.

# Arguments
- `log10pc0`: The log10 first power color. Scalar or array.
- `log10pc1`: The log10 second power color. Scalar or array.

# Keyword Arguments
- `center`: The center coordinates in log10 scale.
  Default: `(log10(4.51920), log10(0.453724))`.

# Returns
- `hue`: Angle(s) in radians.

# References
Heil et al. 2015, MNRAS, 448, 3348.
"""
function hue_from_logpower_color(
    log10pc0, log10pc1;
    center = (log10(4.51920), log10(0.453724)),
)
    return @. (3 / 4) * π - atan(log10pc1 - center[2], log10pc0 - center[1])
end

_limit_angle_to_360(angle) = mod(angle, 360)
