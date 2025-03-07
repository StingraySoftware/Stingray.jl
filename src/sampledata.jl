using DelimitedFiles
using Printf

struct LightCurve
    time::Vector{Float64}
    counts::Vector{Float64}
    dt::Float64
end

function sample_data()
    """
    Import data from .txt file and return a LightCurve object.

    Returns
    -------
    sample: LightCurve object
        The LightCurve object with the desired time stamps
        and counts.
    """
    lc_file = joinpath(@__DIR__, "datasets", "lc_sample.txt")
    data = readdlm(lc_file)

    # Extract first and second columns to indicate dates and counts respectively
    dates = data[:, 1]
    dt = dates[2] - dates[1]
    counts = data[:, 2]

    # Return LightCurve object
    return LightCurve(dates, counts, dt)
end