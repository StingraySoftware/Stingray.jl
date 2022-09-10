@with_kw mutable struct LightCurve{T1<:AbstractVector,T2<:AbstractMatrix}
    time::T1
    counts::Vector{Int}=Int[]
    err::Vector{Float64}=Float64[]
    input_counts::Bool=true
    gti::T2=reshape(Float64[],0,2)
    err_dist::String=""
    mjdref::Float64=0
    dt::Float64=0
    skip_checks::Bool=false
    low_memory::Bool=false
    mission::String=""
    instr::String=""
    header::String=""
end

function make_lightcurves(toa::AbstractVector{<:Real}, dt::Real, 
                          gti::AbstractMatrix{<:Real}; tseg::Real=-1, 
                          tstart::Real=-1, mjdref::Real=0)
    sort!(toa)
    if tstart == -1
        # if tstart is not set, assume light curve starts with first photon
        tstart = toa[begin]
    end

    # compute the number of bins in the light curve
    # for cases where tseg/dt is not integer.
    # TODO: check that this is always consistent and that we
    # are not throwing away good events.

    if tseg == -1
        tseg = toa[end-1] - tstart
    end

    timebin = tseg ÷ dt
    # If we are missing the next bin by just 1%, let's round up:
    if tseg / dt - timebin >= 0.99
        timebin += 1
    end
    
    tend = tstart + timebin * dt
    good = tstart .<= toa .< tend
    
    histbins = range(tstart, tend, step = dt)
    counts= fit(Histogram, keepat!(toa,good), histbins).weights
    time = @view(histbins[begin:end-1]) .+ 0.5 .* dt

    return LightCurve(time=time, counts=counts, gti=gti, mjdref=mjdref, dt=dt,
                      skip_checks=true, err_dist="poisson")
end

function rebin(lc::LightCurve, dt_new::Real; method::String = "sum")
    if dt_new<lc.dt
        throw(ArgumentError("New time resolution must be larger than old time resolution!"))
    end

    bin_time = Float64[]
    bin_counts = Int[]
    bin_err = Float64[]
    gti_new = Vector{Float64}[]

    for g in eachrow(lc.gti)
        if g[2] - g[1] < dt_new
            continue
        else
            # find start and end of GTI segment in data
            start_ind = searchsortedfirst.(Ref(lc.time), g[1])
            end_ind = searchsortedlast.(Ref(lc.time), g[2])

            t_temp = @view(lc.time[start_ind:end_ind])
            c_temp = @view(lc.counts[start_ind:end_ind])

            if !isempty(lc.err)
                e_temp = @view(lc.err[start_ind:end_ind])
            else
                e_temp = Float64[]
            end

            bin_t, bin_c, bin_e, _ = rebin_data(t_temp, c_temp, dt_new,
                                                      y_err=e_temp, method=method)

            append!(bin_time, bin_t)
            append!(bin_counts, bin_c)
            append!(bin_err, bin_e)
            push!(gti_new, [g[1], g[2]])
        end
    end

    if isempty(gti_new) & !isempty(lc.gti)
        throw(ArgumentError("No valid GTIs after rebin."))
    elseif isempty(gti_new)
        gti_new = lc.gti
    else
        gti_new = mapreduce(permutedims, vcat, gti_new)
    end

    return LightCurve(time=bin_time, counts=bin_counts, err=bin_err,
                      mjdref=lc.mjdref, dt=dt_new, gti=gti_new,
                      skip_checks=true)
end
