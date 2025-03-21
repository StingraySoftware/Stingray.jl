using Plots
using Interpolations
using LinearAlgebra
using Colors
using Statistics
using StatsPlots
using DataFrames
module PowerColors

const DEFAULT_COLOR_CONFIGURATION = Dict(
    "center" => [4.51920, 0.453724],
    "ref_angle" => 3 * pi / 4,
    "state_definitions" => Dict(
        "HSS" => Dict("hue_limits" => [300, 360], "color" => "red"),
        "LHS" => Dict("hue_limits" => [-20, 140], "color" => "blue"),
        "HIMS" => Dict("hue_limits" => [140, 220], "color" => "green"),
        "SIMS" => Dict("hue_limits" => [220, 300], "color" => "yellow"),
    ),
    "rms_spans" => Dict(
        -20 => [0.3, 0.7],
        0 => [0.3, 0.7],
        10 => [0.3, 0.6],
        40 => [0.25, 0.4],
        100 => [0.25, 0.35],
        150 => [0.2, 0.3],
        170 => [0.0, 0.3],
        200 => [0, 0.15],
        370 => [0, 0.15],
    ),
)

export DEFAULT_COLOR_CONFIGURATION
end
DEFAULT_COLOR_CONFIGURATION = Dict(
    "center" => [4.51920, 0.453724],
    "ref_angle" => 3 * pi / 4,
    "state_definitions" => Dict(
        "HSS" => Dict("hue_limits" => [300, 360], "color" => "red"),
        "LHS" => Dict("hue_limits" => [-20, 140], "color" => "blue"),
        "HIMS" => Dict("hue_limits" => [140, 220], "color" => "green"),
        "SIMS" => Dict("hue_limits" => [220, 300], "color" => "yellow"),
    ),
    "rms_spans" => Dict(
        -20 => [0.3, 0.7],
        0 => [0.3, 0.7],
        10 => [0.3, 0.6],
        40 => [0.25, 0.4],
        100 => [0.25, 0.35],
        150 => [0.2, 0.3],
        170 => [0.0, 0.3],
        200 => [0, 0.15],
        370 => [0, 0.15],
    ),
)
using Interpolations: LinearInterpolation, Flat
function _get_rms_span_functions(configuration = DEFAULT_COLOR_CONFIGURATION)
    rms_spans = configuration["rms_spans"]

    # Extract and sort unique x-values
    x = sort(unique(collect(keys(rms_spans))))
    ymin = [rms_spans[x_val][1] for x_val in x]
    ymax = [rms_spans[x_val][2] for x_val in x]

    # Create interpolation functions
    ymin_func = LinearInterpolation(x, ymin, extrapolation_bc = Flat())
    ymax_func = LinearInterpolation(x, ymax, extrapolation_bc = Flat())

    return ymin_func, ymax_func
end
function _create_rms_hue_plot(;
    polar = false,
    plot_spans = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    if polar
        fig = plot(proj = :polar)
        ylims!(fig, (0.0, 0.75))
        yticks!(fig, [0, 0.25, 0.5, 0.75, 1])
        plot!(fig, grid = true)
    else
        fig = plot(
            xlims = (0, 360),
            ylims = (0.0, 0.7),
            xlabel = "Hue",
            ylabel = "Fractional rms",
        )
    end
    if !plot_spans
        return fig
    end
    ymin_func, ymax_func = _get_rms_span_functions(configuration)
    for (state, details) in configuration["state_definitions"]
        color = parse(Colorant, details["color"])
        xmin, xmax = details["hue_limits"]

        xs = collect(range(xmin, stop = xmax, length = 20))

        if polar
            # For polar plots, handle each point individually to avoid broadcasting issues
            xs_rad = [deg2rad(x) for x in xs]
            y_mins = [ymin_func(x) for x in xs]
            y_maxs = [ymax_func(x) for x in xs]

            # Plot the minimum line
            plot!(fig, xs_rad, y_mins, color = color, label = "", alpha = 0.5)
            # Plot the maximum line
            plot!(fig, xs_rad, y_maxs, color = color, label = "", alpha = 0.5)
            # Fill between the lines
            for i = 1:(length(xs)-1)
                x_section = [xs_rad[i], xs_rad[i+1], xs_rad[i+1], xs_rad[i]]
                y_section = [y_mins[i], y_mins[i+1], y_maxs[i+1], y_maxs[i]]
                plot!(
                    fig,
                    x_section,
                    y_section,
                    seriestype = :shape,
                    color = color,
                    fillalpha = 0.1,
                    linewidth = 0,
                    label = "",
                )
            end
        else
            # For Cartesian plots, standard approach works
            ymin_vals = [ymin_func(x) for x in xs]
            ymax_vals = [ymax_func(x) for x in xs]
            # Calculate ribbon values explicitly
            ribbon_vals = [ymax - ymin for (ymax, ymin) in zip(ymax_vals, ymin_vals)]
            plot!(
                fig,
                xs,
                ymin_vals,
                ribbon = ribbon_vals,
                fillalpha = 0.1,
                color = color,
                label = "",
            )
        end
    end
    return fig
end

function _limit_angle_to_360(angle)
    while angle >= 360
        angle -= 360
    end
    while angle < 0
        angle += 360
    end
    return angle
end

function _hue_line_data(center, angle; ref_angle = 3 * pi / 4)
    plot_angle = mod(-angle + ref_angle, 2 * pi)

    m = tan(plot_angle)
    if isinf(m)
        x = fill(center[1], 20)
        y = collect(range(-4, stop = 4, length = 20))
    else
        x = collect(range(0, stop = 4, length = 20) .* sign(cos(plot_angle)) .+ center[1])
        y = collect(center[2] .+ m .* (x .- center[1]))
    end
    return x, y
end

function _trace_states(ax, configuration = DEFAULT_COLOR_CONFIGURATION; kwargs...)
    center = log10.(configuration["center"])
    for (state, details) in configuration["state_definitions"]
        color = details["color"]
        hue0, hue1 = details["hue_limits"]
        hue_mean = (hue0 + hue1) / 2
        hue_angle = mod(-deg2rad(hue_mean) + 3 * pi / 4, 2 * pi)

        radius = 1.4
        txt_x = radius * cos(hue_angle) + center[1]
        txt_y = radius * sin(hue_angle) + center[2]
        annotate!(ax, txt_x, txt_y, text(state, :black, :center))

        x0, y0 =
            _hue_line_data(center, deg2rad(hue0); ref_angle = configuration["ref_angle"])
        next_angle = hue0 + 5.0
        x1, y1 =
            _hue_line_data(center, deg2rad(hue0); ref_angle = configuration["ref_angle"])

        while next_angle <= hue1
            x0, y0 = x1, y1
            x1, y1 = _hue_line_data(
                center,
                deg2rad(next_angle);
                ref_angle = configuration["ref_angle"],
            )

            plot!(
                ax,
                [x0[1], x0[end], x1[end], x0[1]],
                [y0[1], y0[end], y1[end], y0[1]],
                seriestype = :shape,
                color = color,
                lw = 0,
                label = "",
            )
            next_angle += 5.0
        end
    end
end

function _create_pc_plot(;
    xrange = (-2, 2),
    yrange = (-2, 2),
    plot_spans = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    fig = plot()
    xlabel!("log_{10}PC1")
    ylabel!("log_{10}PC2")

    if !plot_spans
        xlims!(xrange...)
        ylims!(yrange...)
        return fig
    end

    center = log10.(configuration["center"])
    xlims!(center[1] + xrange[1], center[1] + xrange[2])
    ylims!(center[2] + yrange[1], center[2] + yrange[2])

    for angle = 0:20:360
        x, y =
            _hue_line_data(center, deg2rad(angle); ref_angle = configuration["ref_angle"])
        plot!(x, y, lw = 0.2, ls = :dot, color = :black, alpha = 0.3)
    end

    scatter!([center[1]], [center[2]], marker = :+, color = :black)

    limit_angles = Set(
        vcat(
            [
                configuration["state_definitions"][state]["hue_limits"] for
                state in keys(configuration["state_definitions"])
            ]...,
        ),
    )
    limit_angles = [_limit_angle_to_360(angle) for angle in limit_angles]

    for angle in limit_angles
        x, y =
            _hue_line_data(center, deg2rad(angle); ref_angle = configuration["ref_angle"])
        plot!(x, y, lw = 1, ls = :dot, color = :black, alpha = 1)
    end

    _trace_states(fig, configuration = configuration, alpha = 0.1)

    return fig
end


function plot_power_colors(
    p1::Number,
    p1e::Number,
    p2::Number,
    p2e::Number;
    plot_spans::Bool = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    """
    Plot power colors.

    Parameters
    ----------
    p1 : Number
        The first power value.
    p1e : Number
        The error in the first power value.
    p2 : Number
        The second power value.
    p2e : Number
        The error in the second power value.

    Other Parameters
    ----------------
    center : (Number, Number)
        The center coordinates of the plot. Default is (4.51920, 0.453724).
    plot_spans : Bool
        Whether to plot the spans. Default is false.
    configuration: Dict
        The color configuration to use. Default is DEFAULT_COLOR_CONFIGURATION.

    Returns
    -------
    fig : Plot
        The Julia Plot object representing the power color plot.
    """
    p1e = 1 / p1 * p1e
    p2e = 1 / p2 * p2e
    p1 = log10(p1)
    p2 = log10(p2)

    fig = _create_pc_plot(plot_spans = plot_spans, configuration = configuration)

    scatter!([p1], [p2], marker = :circle, color = :black)

    df = DataFrame(x = [p1], y = [p2], xerr = [p1e], yerr = [p2e])

    # Use error bars correctly without `permute`
    plot!(df.x, df.y, ribbon = (df.yerr, df.yerr), lw = 0, alpha = 0.4, color = :black)
    plot!(df.x, df.y, ribbon = (df.xerr, df.xerr), lw = 0, alpha = 0.4, color = :black)

    return fig
end

function hue_from_power_color(pc1, pc2; center = nothing)
    if center === nothing
        return atan.(pc2, pc1)
    else
        # Convert to log space
        log_pc1 = log10(pc1)
        log_pc2 = log10(pc2)
        log_center = log10.(center)

        # Calculate the raw angle from the center in log space
        raw_angle = atan(log_pc2 - log_center[2], log_pc1 - log_center[1])
        return 3 / 4 * π - raw_angle
    end
end
function plot_hues(
    rms,
    rmse,
    pc1,
    pc2;
    polar = false,
    plot_spans = false,
    configuration = DEFAULT_COLOR_CONFIGURATION,
)
    """
    Plot hues.

    Parameters
    ----------
    rms : Array{Number}
        The root mean square values.
    rmse : Array{Number}
        The errors in the root mean square values.
    pc1 : Array{Number}
        The first power color component.
    pc2 : Array{Number}
        The second power color component.

    Other Parameters
    ----------------
    polar : Bool
        Whether to use a polar plot. Default is false.
    plot_spans : Bool
        Whether to plot the spans. Default is false.
    configuration: Dict
        The color configuration to use. Default is DEFAULT_COLOR_CONFIGURATION.

    Returns
    -------
    ax : Plot
        The Julia Plot object representing the hues plot.
    """
    hues = hue_from_power_color(collect(pc1), collect(pc2))  # Ensure arrays

    ax = _create_rms_hue_plot(
        polar = polar,
        plot_spans = plot_spans,
        configuration = configuration,
    )
    hues = mod.(hues, 2 * pi)
    if !polar
        hues = rad2deg.(hues)
    end

    scatter!(hues, rms, yerror = rmse, alpha = 0.5)  # Correct error bars

    return ax
end

function integrate_power_in_frequency_range(
    frequency,
    power,
    freq_range;
    power_err = nothing,
    df = nothing,
    m = 1,
    poisson_power = 0,
)
    """
    Integrate the power over a frequency range.

    Parameters
    ----------
    frequency : Vector{Number}
        The frequency array.
    power : Vector{Number}
        The power at each frequency.
    freq_range : Tuple{Number, Number} or Vector{Number}
        The frequency range over which to integrate.

    Other Parameters
    ----------------
    power_err : Vector{Number}, optional
        The error on the power.
    df : Number, optional
        The frequency resolution.
    m : Int, optional
        The number of power spectrum averages.
    poisson_power : Number, optional
        The Poisson noise level.

    Returns
    -------
    variance : Number
        The integrated power (variance) in the frequency range.
    variance_err : Number
        The error on the integrated power.
    """
    frequency = collect(frequency)
    power = collect(power)

    # Ensure freq_range is a two-element vector
    if length(freq_range) != 2
        throw(ArgumentError("freq_range must have exactly 2 elements"))
    end

    f0, f1 = freq_range

    # Calculate frequency resolution if not provided
    df = isnothing(df) ? median(diff(frequency)) : df

    # Define frequency range including bin width
    input_frequency_low_edges = frequency .- df / 2
    input_frequency_high_edges = frequency .+ df / 2

    # Find indices of frequencies within the required range
    idx = findall((input_frequency_high_edges .>= f0) .& (input_frequency_low_edges .<= f1))

    if isempty(idx)
        throw(ArgumentError("No frequencies found in the range ($f0, $f1)"))
    end

    # Get power values within range
    power_in_range = power[idx]

    # Handle errors if provided
    power_err_in_range = isnothing(power_err) ? power_in_range ./ sqrt(m) : power_err[idx]

    # Calculate total variance (integrated power)
    # For power spectral density, integration is multiplication by df
    variance = sum(power_in_range .- poisson_power) * df

    # Error propagation for the sum
    variance_err = sqrt(sum(power_err_in_range .^ 2)) * df

    return variance, variance_err
end

function power_color(
    frequency,
    power;
    power_err = nothing,
    freq_edges = [1 / 256, 1 / 32, 0.25, 2.0, 16.0],
    df = nothing,
    m = 1,
    freqs_to_exclude = nothing,
    poisson_power = 0,
    return_log = false,
)
    """
    Calculate two power colors from a power spectrum.

    Power colors are an alternative to spectral colors to understand the spectral state of an
    accreting source. They are defined as the ratio of the power in two frequency ranges,
    analogously to the colors calculated from electromagnetic spectra.

    This function calculates two power colors, using the four frequency ranges contained
    between the five frequency edges in `freq_edges`. Given [f0, f1, f2, f3, f4], the
    two power colors are calculated as the following ratios of the integrated power
    (which are variances):

    - PC0 = Var([f0, f1]) / Var([f2, f3])
    - PC1 = Var([f1, f2]) / Var([f3, f4])

    Errors are calculated using simple error propagation from the integrated power errors.

    See Heil et al. 2015, MNRAS, 448, 3348.

    Parameters
    ----------
    frequency : Vector{Number}
        The frequencies of the power spectrum.
    power : Vector{Number}
        The power at each frequency.

    Other Parameters
    ----------------
    power_err : Vector{Number}, optional
        The power error bar at each frequency.
    freq_edges : Vector{Number}, optional
        The five edges defining the four frequency intervals to use to calculate the power color.
    df : Number, optional
        The frequency resolution of the input data. If `nothing`, it is calculated from the median difference of input frequencies.
    m : Int, optional
        The number of segments and/or contiguous frequency bins averaged to obtain power.
    freqs_to_exclude : Vector{Tuple{Number, Number}}, optional
        The ranges of frequencies to exclude from the calculation of the power color.
    poisson_power : Number, optional
        The Poisson noise level of the power spectrum.
    return_log : Bool, optional
        Return the base-10 logarithm of the power color and the errors.

    Returns
    -------
    PC0 : Number
        The first power color.
    PC0_err : Number
        The error on the first power color.
    PC1 : Number
        The second power color.
    PC1_err : Number
        The error on the second power color.
    """
    freq_edges = collect(freq_edges)
    length(freq_edges) == 5 || throw(ArgumentError("freq_edges must have 5 elements"))

    frequency = collect(frequency)
    power = collect(power)

    df = isnothing(df) ? median(diff(frequency)) : df
    input_frequency_low_edges = frequency .- df / 2
    input_frequency_high_edges = frequency .+ df / 2

    # Fix the error check - ensure freq_edges are within the range of input frequencies
    if minimum(freq_edges) < minimum(input_frequency_low_edges)
        throw(
            ArgumentError(
                "The minimum frequency edge is lower than the available frequencies",
            ),
        )
    end

    if maximum(freq_edges) > maximum(input_frequency_high_edges)
        throw(
            ArgumentError(
                "The maximum frequency edge is higher than the available frequencies",
            ),
        )
    end

    power_err = isnothing(power_err) ? power ./ sqrt(m) : collect(power_err)

    if !isnothing(freqs_to_exclude)
        for (f0, f1) in freqs_to_exclude
            frequency_mask =
                (input_frequency_low_edges .> f0) .& (input_frequency_high_edges .< f1)
            idx0, idx1 = searchsortedfirst(frequency, f0), searchsortedlast(frequency, f1)
            power[frequency_mask] .= mean([power[idx0], power[idx1]])
        end
    end

    var00, var00_err = integrate_power_in_frequency_range(
        frequency,
        power,
        freq_edges[1:2],
        power_err = power_err,
        df = df,
        m = m,
        poisson_power = poisson_power,
    )
    var01, var01_err = integrate_power_in_frequency_range(
        frequency,
        power,
        freq_edges[3:4],
        power_err = power_err,
        df = df,
        m = m,
        poisson_power = poisson_power,
    )
    var10, var10_err = integrate_power_in_frequency_range(
        frequency,
        power,
        freq_edges[2:3],
        power_err = power_err,
        df = df,
        m = m,
        poisson_power = poisson_power,
    )
    var11, var11_err = integrate_power_in_frequency_range(
        frequency,
        power,
        freq_edges[4:5],
        power_err = power_err,
        df = df,
        m = m,
        poisson_power = poisson_power,
    )

    pc0 = var00 / var01
    pc1 = var10 / var11
    pc0_err = pc0 * (var00_err / var00 + var01_err / var01)
    pc1_err = pc1 * (var10_err / var10 + var11_err / var11)

    if return_log
        pc0_err = 1 / pc0 * pc0_err
        pc1_err = 1 / pc1 * pc1_err
        pc0 = log10(pc0)
        pc1 = log10(pc1)
    end

    return pc0, pc0_err, pc1, pc1_err
end

function hue_from_logpower_color(log10pc0, log10pc1; center = log10.([4.51920, 0.453724]))
    """
    Measure the angle of a point in the log-power color diagram with respect to the center.

    Angles are measured in radians, **in the clockwise direction**, with respect to a line oriented
    at -45 degrees with respect to the horizontal axis.

    See Heil et al. 2015, MNRAS, 448, 3348.

    Parameters
    ----------
    log10pc0 : Number
        The log10 power color in the first frequency range.
    log10pc1 : Number
        The log10 power color in the second frequency range.

    Other Parameters
    ----------------
    center : Vector{Number}, optional, default log10([4.51920, 0.453724])
        The coordinates of the center of the power color diagram.

    Returns
    -------
    hue : Number
        The angle of the point with respect to the center, in radians.
    """
    return (3 / 4) * pi - atan(log10pc1 - center[2], log10pc0 - center[1])
end
function polar_rms_scatter_with_hues(
    state_definitions,
    pc1_vals,
    pc2_vals,
    rms_vals;
    title = "Polar Power Color Hue Plot",
)
    """
    Generate a polar scatter plot of power color hues with background state spans.

    Parameters
    ----------
    state_definitions : Dict
        Dictionary containing hue limits and colors for different states.
    pc1_vals : Vector{Float64}
        First power color component values.
    pc2_vals : Vector{Float64}
        Second power color component values.
    rms_vals : Vector{Float64}
        Corresponding fractional RMS values.

    Returns
    -------
    fig : Plot
        The generated polar scatter plot with hues.
    """

    # Compute hues from power color components
    hues_rad = atan.(pc2_vals, pc1_vals)  # Compute hue as the angle from power color components
    hues_deg = rad2deg.(hues_rad)  # Convert to degrees

    # Ensure hues are in the range [0, 360]
    hues_deg = mod.(hues_deg, 360.0)
    hues_rad = deg2rad.(hues_deg)

    # Create Polar Plot
    fig = plot(proj = :polar, legend = false, title = title, ylim = (0, 0.75))

    # Plot background state spans
    for (state, params) in state_definitions
        hue_min, hue_max = params["hue_limits"]
        color = params["color"]

        # Convert hue limits to radians
        theta1, theta2 = deg2rad(hue_min), deg2rad(hue_max)

        # Define radial (fractional RMS) span
        r_vals = [0.0, 0.75, 0.75, 0.0]
        theta_vals = [theta1, theta1, theta2, theta2]

        # Plot the background span
        plot!(theta_vals, r_vals, fillalpha = 0.15, lw = 0, color = color, label = "")
    end

    # Scatter plot for hues (Power Color representation)
    scatter!(
        hues_rad,
        rms_vals,
        markersize = 6,
        color = :blue,
        alpha = 0.8,
        markerstrokewidth = 0,
    )

    return fig
end
