using Plots
using Interpolations
using LinearAlgebra
using colors

DEFAULT_COLOR_CONFIGURATION = Dict(
    "center" => [4.51920, 0.453724],
    "ref_angle" => 3 * pi / 4,
    "state_definitions" => Dict(
        "HSS" => Dict("hue_limits" => [300, 360], "color" => "red"),
        "LHS" => Dict("hue_limits" => [-20, 140], "color" => "blue"),
        "HIMS" => Dict("hue_limits" => [140, 220], "color" => "green"),
        "SIMS" => Dict("hue_limits" => [220, 300], "color" => "yellow")
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
        370 => [0, 0.15]
    )
)

using Interpolations

function _get_rms_span_functions(configuration=DEFAULT_COLOR_CONFIGURATION)
    rms_spans = configuration["rms_spans"]

    x = collect(keys(rms_spans))
    ymin = [v[1] for v in values(rms_spans)]
    ymax = [v[2] for v in values(rms_spans)]

    ymin_func = LinearInterpolation(x, ymin)
    ymax_func = LinearInterpolation(x, ymax)

    return ymin_func, ymax_func
end

function _create_rms_hue_plot(; polar=false, plot_spans=false, configuration=DEFAULT_COLOR_CONFIGURATION)
    if polar
        fig = plot(proj=:polar)
        rlims!(fig, (0, 0.75))
        yticks!(fig, [0, 0.25, 0.5, 0.75, 1])
        grid!(fig, true)
    else
        fig = plot(xlims=(0, 360), ylims=(0, 0.7), xlabel="Hue", ylabel="Fractional rms")
    end

    if !plot_spans
        return fig
    end

    ymin_func, ymax_func = _get_rms_span_functions(configuration)

    for (state, details) in configuration["state_definitions"]
        color = parse(Colorant, details["color"])
        xmin, xmax = details["hue_limits"]

        x_func = x -> x
        if polar
            x_func = x -> deg2rad(x)
        end

        xs = range(x_func(xmin), stop=x_func(xmax), length=20)
        ymin_vals = ymin_func(range(xmin, stop=xmax, length=20))
        ymax_vals = ymax_func(range(xmin, stop=xmax, length=20))

        plot!(fig, xs, ymin_vals, fill_between=(ymin_vals, ymax_vals), fillalpha=0.1, color=color, label="")

        if xmin < 0 && !polar
            xs_extra = range(x_func(xmin + 360), stop=x_func(360), length=20)
            ymin_vals_extra = ymin_func(range(xmin, stop=0, length=20))
            ymax_vals_extra = ymax_func(range(xmin, stop=0, length=20))
            plot!(fig, xs_extra, ymin_vals_extra, fill_between=(ymin_vals_extra, ymax_vals_extra), fillalpha=0.1, color=color, label="")
        end
        if xmax > 360 && !polar
            xs_extra = range(x_func(0), stop=x_func(xmax - 360), length=20)
            ymin_vals_extra = ymin_func(range(360, stop=xmax, length=20))
            ymax_vals_extra = ymax_func(range(360, stop=xmax, length=20))
            plot!(fig, xs_extra, ymin_vals_extra, fill_between=(ymin_vals_extra, ymax_vals_extra), fillalpha=0.1, color=color, label="")
        end
    end

    return fig
end
