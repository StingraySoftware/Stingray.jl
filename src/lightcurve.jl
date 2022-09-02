@with_kw mutable struct LightCurve{T1<:AbstractVector,T2<:AbstractMatrix}
    time::T1
    counts::Int=0
    err::Vector{Float64}=Float64[]
    input_counts::Bool=true
    gti::T2=reshape(Float64[],0,2)
    err_dist::String=""
    mjdref::Float64=0
    skip_checks::Bool=false
    low_memory::Bool=false
    mission::String=""
    instr::String=""
    header::String=""
end

function make_lightcurves(toa::AbstractVector{<:Real}, dt::Float64; tseg::Real=0.0,
                          tstart::Real=nothing, gti::AbstractMatrix{<:Real}=nothing,
                          mjdref::Real=0, use_hist::Bool=false)
    sort(toa)
    if isnothing(tstart)
        # if tstart is not set, assume light curve starts with first photon
        tstart = toa[begin]
    end

    # compute the number of bins in the light curve
    # for cases where tseg/dt is not integer.
    # TODO: check that this is always consistent and that we
    # are not throwing away good events.

    if isnothing(tseg)
        tseg = toa[end-1] - tstart
    end

    timebin = tseg ÷ dt
    # If we are missing the next bin by just 1%, let's round up:
    if tseg / dt - timebin >= 0.99
        timebin += 1
    end
    
    tend = tstart + timebin * dt
    good = tstart <= toa < tend
    if !use_hist
        binned_toas = (keepat!(toa,good) .- tstart) .÷ dt
        counts = bincount(binned_toas, minlength=timebin)
        time = tstart + range(0.5, 0.5 + len(counts)) * dt
    else
        histbins = range(tstart, tend + dt, dt)
        counts, histbins = histogram(keepat!(toa,good), bins=histbins)
        time = @view(histbins[begin:end-1]) + 0.5 * dt
    end

    return Lightcurve(time, counts; gti=gti, mjdref=mjdref, dt=dt,
                      skip_checks=True, err_dist="poisson")
end