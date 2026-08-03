# ──────────────────────────────────────────────────────────────────────────────
# Lomb-Scargle Periodograms for unevenly sampled data
#
# Mirrors Python Stingray's stingray/lombscargle.py
# ──────────────────────────────────────────────────────────────────────────────

"""
    autofrequency(; min_freq=nothing, max_freq=nothing, df=nothing,
                    dt=nothing, length=nothing, nyquist_factor=1)

Build the frequency grid for a Lomb-Scargle periodogram.

Determines `df`, `max_freq`, and `min_freq` from whatever subset of
parameters the user provides, then returns a regularly spaced frequency
array from `min_freq` to `max_freq` in steps of `df`.

Mirrors Python Stingray's `lombscargle._autofrequency`.

# Keyword Arguments
- `min_freq::Union{Nothing, Real}=nothing`: Minimum frequency. Defaults to `df / 2`.
- `max_freq::Union{Nothing, Real}=nothing`: Maximum frequency. Defaults to Nyquist: `nyquist_factor * 0.5 / dt`.
- `df::Union{Nothing, Real}=nothing`: Frequency resolution. Defaults to `1 / length`.
- `dt::Union{Nothing, Real}=nothing`: Time resolution of the data.
- `length::Union{Nothing, Real}=nothing`: Total duration of the observation.
- `nyquist_factor::Real=1`: Multiplier on the Nyquist frequency.

# Returns
- `Vector{Float64}`: The array of frequencies.

# Examples
```julia
freqs = autofrequency(min_freq=0.1, max_freq=0.5, df=0.1)
# [0.1, 0.2, 0.3, 0.4, 0.5]

freqs = autofrequency(min_freq=0.1, dt=1.0, length=10.0)
# [0.1, 0.2, 0.3, 0.4, 0.5]
```
"""
function autofrequency(; min_freq::Union{Nothing, Real}=nothing,
                         max_freq::Union{Nothing, Real}=nothing,
                         df::Union{Nothing, Real}=nothing,
                         dt::Union{Nothing, Real}=nothing,
                         length::Union{Nothing, Real}=nothing,
                         nyquist_factor::Real=1)

    # Determine df
    if (isnothing(df) || df <= 0) && isnothing(length)
        throw(ValueError("Either df or length must be specified."))
    elseif isnothing(df) || df <= 0
        df = 1.0 / length
    end

    # Determine max_freq
    if isnothing(max_freq) && (isnothing(dt) || dt == 0)
        throw(ValueError("Either max_freq or dt must be specified."))
    elseif isnothing(max_freq)
        max_freq = nyquist_factor * 0.5 / dt
    end

    # Determine min_freq
    if isnothing(min_freq)
        min_freq = df / 2.0
    elseif min_freq <= 0
        @warn "min_freq must be positive and >0. Setting to df / 2."
        min_freq = df / 2.0
    end

    freq = collect(Float64(min_freq):Float64(df):Float64(max_freq + df))
    # Trim any values beyond max_freq (floating point overshoot)
    filter!(f -> f <= max_freq + df * 0.5, freq)

    return freq
end

# ValueError alias for consistency with Python-style error messages
struct ValueError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ValueError) = print(io, "ValueError: ", e.msg)
