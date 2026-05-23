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

function _get_rms_span_functions(configuration = DEFAULT_COLOR_CONFIGURATION)
    rms_spans = configuration["rms_spans"]
    x = sort(collect(keys(rms_spans)))
    ymin = [rms_spans[k][1] for k in x]
    ymax = [rms_spans[k][2] for k in x]

    x_float = Float64.(x)
    ymin_func = linear_interpolation(x_float, Float64.(ymin))
    ymax_func = linear_interpolation(x_float, Float64.(ymax))
    return ymin_func, ymax_func
end

function _hue_line_data(center, angle; ref_angle = 3 * π / 4)
    plot_angle = mod(-angle + ref_angle, 2π)
    m = tan(plot_angle)

    if isinf(m)
        x = fill(center[1], 20)
        y = range(-4, 4; length = 20)
    else
        x = range(0, 4; length = 20) .* sign(cos(plot_angle)) .+ center[1]
        y = center[2] .+ m .* (x .- center[1])
    end
    return collect(x), collect(y)
end

function _create_rms_hue_plot(;
    polar::Bool = false,
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    fig = CairoMakie.Figure(size = (600, 500))

    if polar
        ax = CairoMakie.PolarAxis(fig[1, 1]; rlimits = (0, 0.75),
            rticks = [0, 0.25, 0.5, 0.75])
    else
        ax = CairoMakie.Axis(fig[1, 1];
            xlabel = "Hue (degrees)", ylabel = "Fractional rms",
            limits = (0, 360, 0, 0.7))
    end

    if !plot_spans
        return fig, ax
    end

    ymin_func, ymax_func = _get_rms_span_functions(configuration)

    for (state, params) in configuration["state_definitions"]
        color = params["color"]
        xmin, xmax = params["hue_limits"]

        x_lin = range(xmin, xmax; length = 20)
        y_low = [ymin_func(xi) for xi in x_lin]
        y_high = [ymax_func(xi) for xi in x_lin]

        x_plot = polar ? deg2rad.(collect(x_lin)) : collect(x_lin)

        CairoMakie.band!(ax, x_plot, y_low, y_high; color = (Symbol(color), 0.1))

        if !polar && xmin < 0
            x_wrap = range(xmin + 360, 360; length = 20)
            x_orig = range(xmin, 0; length = 20)
            y_low_w = [ymin_func(xi) for xi in x_orig]
            y_high_w = [ymax_func(xi) for xi in x_orig]
            CairoMakie.band!(ax, collect(x_wrap), y_low_w, y_high_w;
                color = (Symbol(color), 0.1))
        end
        if !polar && xmax > 360
            x_wrap = range(0, xmax - 360; length = 20)
            x_orig = range(360, xmax; length = 20)
            y_low_w = [ymin_func(xi) for xi in x_orig]
            y_high_w = [ymax_func(xi) for xi in x_orig]
            CairoMakie.band!(ax, collect(x_wrap), y_low_w, y_high_w;
                color = (Symbol(color), 0.1))
        end
    end

    return fig, ax
end

function _trace_states(ax, configuration = DEFAULT_COLOR_CONFIGURATION; alpha = 0.1)
    center = log10.(configuration["center"])

    for (state, params) in configuration["state_definitions"]
        color = params["color"]
        hue0, hue1 = params["hue_limits"]
        hue_mean = (hue0 + hue1) / 2
        hue_angle = mod(-deg2rad(hue_mean) + 3π / 4, 2π)

        radius = 1.4
        txt_x = radius * cos(hue_angle) + center[1]
        txt_y = radius * sin(hue_angle) + center[2]
        CairoMakie.text!(ax, txt_x, txt_y; text = state, align = (:center, :center),
            color = :black, fontsize = 12)

        next_angle = hue0 + 5.0
        x0, y0 = _hue_line_data(center, deg2rad(hue0);
            ref_angle = configuration["ref_angle"])

        while next_angle <= hue1
            x1, y1 = _hue_line_data(center, deg2rad(next_angle);
                ref_angle = configuration["ref_angle"])
            tri_x = [x0[1], x0[end], x1[end], x0[1]]
            tri_y = [y0[1], y0[end], y1[end], y0[1]]
            CairoMakie.poly!(ax, CairoMakie.Point2f.(zip(tri_x, tri_y));
                color = (Symbol(color), alpha), strokewidth = 0)
            x0, y0 = x1, y1
            next_angle += 5.0
        end
    end
end

function _create_pc_plot(;
    xrange = [-2, 2],
    yrange = [-2, 2],
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    fig = CairoMakie.Figure(size = (600, 600))

    if !plot_spans
        ax = CairoMakie.Axis(fig[1, 1];
            xlabel = "log₁₀PC1", ylabel = "log₁₀PC2", aspect = 1,
            limits = (xrange[1], xrange[2], yrange[1], yrange[2]))
        return fig, ax
    end

    center = log10.(configuration["center"])
    ax = CairoMakie.Axis(fig[1, 1];
        xlabel = "log₁₀PC1", ylabel = "log₁₀PC2", aspect = 1,
        limits = (center[1] + xrange[1], center[1] + xrange[2],
                  center[2] + yrange[1], center[2] + yrange[2]))

    for angle in 0:20:359
        x, y = _hue_line_data(center, deg2rad(angle);
            ref_angle = configuration["ref_angle"])
        CairoMakie.lines!(ax, x, y; linewidth = 0.2, linestyle = :dot,
            color = (:black, 0.3))
    end

    CairoMakie.scatter!(ax, [center[1]], [center[2]]; marker = '+',
        markersize = 15, color = :black)

    limit_angles = Set{Float64}()
    for (_, params) in configuration["state_definitions"]
        for a in params["hue_limits"]
            push!(limit_angles, _limit_angle_to_360(Float64(a)))
        end
    end

    for angle in limit_angles
        x, y = _hue_line_data(center, deg2rad(angle);
            ref_angle = configuration["ref_angle"])
        CairoMakie.lines!(ax, x, y; linewidth = 1, linestyle = :dot,
            color = (:black, 1.0))
    end

    _trace_states(ax, configuration; alpha = 0.1)

    return fig, ax
end

"""
    plot_power_colors(p1, p1e, p2, p2e; plot_spans=false, configuration=DEFAULT_COLOR_CONFIGURATION)

Plot power colors in the log₁₀PC1 vs log₁₀PC2 plane using CairoMakie.

# Returns
`(fig, ax)` — CairoMakie Figure and Axis objects.
"""
function plot_power_colors(
    p1, p1e, p2, p2e;
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    p1e_log = (1 / p1) * p1e
    p2e_log = (1 / p2) * p2e
    p1_log = log10(p1)
    p2_log = log10(p2)

    fig, ax = _create_pc_plot(; plot_spans = plot_spans, configuration = configuration)
    CairoMakie.errorbars!(ax, [p1_log], [p2_log], [p1e_log]; direction = :x,
        color = (:black, 0.4))
    CairoMakie.errorbars!(ax, [p1_log], [p2_log], [p2e_log]; direction = :y,
        color = (:black, 0.4))
    CairoMakie.scatter!(ax, [p1_log], [p2_log]; color = :black, markersize = 8)

    return fig, ax
end

"""
    plot_hues(rms, rmse, pc1, pc2; polar=false, plot_spans=false, configuration=DEFAULT_COLOR_CONFIGURATION)

Plot hue angle vs fractional rms using CairoMakie.

# Returns
`(fig, ax)` — CairoMakie Figure and Axis objects.
"""
function plot_hues(
    rms, rmse, pc1, pc2;
    polar::Bool = false,
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    pc1_arr = pc1 isa Real ? [pc1] : collect(pc1)
    pc2_arr = pc2 isa Real ? [pc2] : collect(pc2)
    rms_arr = rms isa Real ? [rms] : collect(rms)
    rmse_arr = rmse isa Real ? [rmse] : collect(rmse)

    hues = hue_from_power_color(pc1_arr, pc2_arr)
    hues = mod.(hues, 2π)

    fig, ax = _create_rms_hue_plot(; polar = polar, plot_spans = plot_spans,
        configuration = configuration)

    hue_plot = polar ? hues : rad2deg.(hues)

    CairoMakie.errorbars!(ax, hue_plot, rms_arr, rmse_arr; color = (:black, 0.5))
    CairoMakie.scatter!(ax, hue_plot, rms_arr; markersize = 8, color = (:steelblue, 0.8))

    return fig, ax
end

# RecipesBase / Plots.jl support

"""
    PowerColorPlot(pc0, pc0_err, pc1, pc1_err; plot_spans=false, configuration=DEFAULT_COLOR_CONFIGURATION)

Wrapper type for plotting power colors via `Plots.jl`.

# Usage
```julia
using Plots
plot(PowerColorPlot(pc0, pc0_err, pc1, pc1_err))
plot(PowerColorPlot(pc0, pc0_err, pc1, pc1_err; plot_spans=true))
```
"""
struct PowerColorPlot
    pc0
    pc0_err
    pc1
    pc1_err
    plot_spans::Bool
    configuration::Dict{String, Any}
end

function PowerColorPlot(pc0, pc0_err, pc1, pc1_err;
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    return PowerColorPlot(pc0, pc0_err, pc1, pc1_err, plot_spans, configuration)
end

@recipe function f(p::PowerColorPlot)
    pc0_log = log10(p.pc0)
    pc1_log = log10(p.pc1)
    pc0_err_log = (1 / p.pc0) * p.pc0_err
    pc1_err_log = (1 / p.pc1) * p.pc1_err

    xlabel --> "log₁₀PC1"
    ylabel --> "log₁₀PC2"
    aspect_ratio --> :equal
    legend --> false

    if p.plot_spans
        center = log10.(p.configuration["center"])
        xlims --> (center[1] - 2, center[1] + 2)
        ylims --> (center[2] - 2, center[2] + 2)
    end

    @series begin
        seriestype := :scatter
        markersize --> 5
        markercolor --> :black
        xerror := [pc0_err_log]
        yerror := [pc1_err_log]
        [pc0_log], [pc1_log]
    end
end

"""
    HuePlot(rms, rms_err, pc0, pc1; polar=false, plot_spans=false, configuration=DEFAULT_COLOR_CONFIGURATION)

Wrapper type for plotting hue vs rms via `Plots.jl`.

# Usage
```julia
using Plots
plot(HuePlot(rms, rms_err, pc0, pc1))
plot(HuePlot(rms, rms_err, pc0, pc1; plot_spans=true))
```
"""
struct HuePlot
    rms
    rms_err
    pc0
    pc1
    polar::Bool
    plot_spans::Bool
    configuration::Dict{String, Any}
end

function HuePlot(rms, rms_err, pc0, pc1;
    polar::Bool = false,
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    return HuePlot(rms, rms_err, pc0, pc1, polar, plot_spans, configuration)
end

@recipe function f(h::HuePlot)
    pc0_arr = h.pc0 isa Real ? [h.pc0] : collect(h.pc0)
    pc1_arr = h.pc1 isa Real ? [h.pc1] : collect(h.pc1)
    rms_arr = h.rms isa Real ? [h.rms] : collect(h.rms)
    rmse_arr = h.rms_err isa Real ? [h.rms_err] : collect(h.rms_err)

    hues = hue_from_power_color(pc0_arr, pc1_arr)
    hues = mod.(hues, 2π)
    hue_plot = h.polar ? hues : rad2deg.(hues)

    xlabel --> (h.polar ? "Hue (rad)" : "Hue (degrees)")
    ylabel --> "Fractional rms"
    legend --> false

    if !h.polar
        xlims --> (0, 360)
        ylims --> (0, 0.7)
    end

    @series begin
        seriestype := :scatter
        markersize --> 5
        markercolor --> :steelblue
        markeralpha --> 0.8
        yerror := rmse_arr
        hue_plot, rms_arr
    end
end
