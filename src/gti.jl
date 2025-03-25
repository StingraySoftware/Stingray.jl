function get_total_gti_length(gti::AbstractMatrix{<:Real}; minlen::Real=0.0)
    lengths = diff(gti; dims =2)
    return sum(x->x > minlen ? x : zero(x), lengths)
end

function load_gtis(fits_file::String, gtistring::String="GTI")
    gti = FITS(fits_file) do lchdulist
        gtihdu = lchdulist[gtistring]
        get_gti_from_hdu(gtihdu)
    end
    return gti
end

function get_gti_from_hdu(gtihdu::TableHDU)

    if "START" in FITSIO.colnames(gtihdu)
        startstr = "START"
        stopstr = "STOP"
    else
        startstr = "Start"
        stopstr = "Stop"
    end

    gtistart = read(gtihdu,startstr)
    gtistop = read(gtihdu,stopstr)

    return mapreduce(permutedims, vcat, 
    [[a, b] for (a,b) in zip(gtistart, gtistop)])
end

function check_gtis(gti::AbstractMatrix)

    if ndims(gti) != 2 || size(gti,2) != 2
        throw(ArgumentError("Please check the formatting of the GTIs. 
       They need to be provided as [[gti00 gti01]; [gti10 gti11]; ...]."))
    end

    gti_start = @view gti[:, 1]
    gti_end = @view gti[:, 2]

    if any(gti_end < gti_start)
        throw(ArgumentError(
            "The GTI end times must be larger than the GTI start times."
        )) 
    end

    if any(@view(gti_start[begin+1:end]) < @view(gti_end[begin:end-1]))
        throw(ArgumentError(
            "This GTI has overlaps"
        ))
    end
end

function create_gti_mask(times::AbstractVector{<:Real},gtis::AbstractMatrix{<:Real};
                         safe_interval::AbstractVector{<:Real}=[0,0], min_length::Real=0,
                         dt::Real = -1, epsilon::Real = 0.001)

    if isempty(times)
        throw(ArgumentError("Passing an empty time array to create_gti_mask"))
    end

    check_gtis(gtis)
    mask = zeros(Bool,length(times))

    if min_length>0
        gtis = gtis[min_length .< @view(gtis[:,2]) - @view(gtis[:,1]),:]
            
        if size(gtis,1) < 1
            @warn "No GTIs longer than min_length $(min_length)"
            return mask, gtis
        end
    end   

    if dt < 0
        dt = Statistics.median(diff(times))
    end
    epsilon_times_dt = epsilon * dt

    new_gtis = [[0.0, 0.0] for _ in range(1,size(gtis,1))]
    new_gti_mask = zeros(Bool, size(gtis,1))

    gti_start = @view gtis[:, 1]
    gti_end = @view gtis[:, 2]

    for (ig,(limmin,limmax)) in enumerate(zip(gti_start,gti_end))
        limmin += safe_interval[1]
        limmax -= safe_interval[2]
        if limmax - limmin >= min_length
            new_gtis[ig][:] .= limmin, limmax
            for (i,t) in enumerate(times) 
                if (limmin + dt / 2 - epsilon_times_dt) <= t <= (limmax - dt / 2 + epsilon_times_dt)
                    mask[i] = true
                end
            end
            new_gti_mask[ig] = true
        end
    end

    return mask, mapreduce(permutedims, vcat, keepat!(new_gtis,new_gti_mask))
end

function create_gti_from_condition(time::AbstractVector{<:Real}, condition::AbstractVector{Bool};
    safe_interval::AbstractVector{<:Real}=[0,0], dt::AbstractVector{<:Real}=Float64[])
    
    if length(time) != length(condition)
        throw(ArgumentError("The length of the condition and time arrays must be the same."))
    end

    idxs = contiguous_regions(condition)

    if isempty(dt)
        dt = zero(time) .+ (time[2] .- time[1]) ./ 2
    end

    gtis = Vector{Float64}[]
    for idx in eachrow(idxs)
        startidx = idx[1]
        stopidx = idx[2] - 1

        t0 = time[startidx] - dt[startidx] + safe_interval[1]
        t1 = time[stopidx] + dt[stopidx] - safe_interval[2]
        if t1 - t0 < 0
            continue
        end
        push!(gtis,[t0, t1])
    end
    return mapreduce(permutedims, vcat, gtis)
end

function operations_on_gtis(gti_list::AbstractVector{<:AbstractMatrix{T}}, 
                            operation::Function) where {T<:Real}

    required_interval = nothing

    for gti in gti_list
        check_gtis(gti)

        combined_gti = Interval{T}[]
        for ig in eachrow(gti)
            push!(combined_gti,Interval{Closed,Open}(ig[1],ig[2]))
        end
        if isnothing(required_interval)
            required_interval = IntervalSet(combined_gti)
        else
            required_interval = operation(required_interval, IntervalSet(combined_gti))
        end
    end

    final_gti = Vector{T}[]

    for interval in required_interval.items
        push!(final_gti, [first(interval), last(interval)])
    end

    return mapreduce(permutedims, vcat, final_gti)
end

function get_btis(gtis::AbstractMatrix{<:Real})
    if isempty(gtis)
        throw(ArgumentError("Empty GTI and no valid start_time and stop_time"))
    end
    return get_btis(gtis, gtis[1,1], gtis[end,2])
end

function get_btis(gtis::AbstractMatrix{T}, start_time, stop_time) where {T<:Real}
    if isempty(gtis)
        return T[start_time stop_time]
    end
    check_gtis(gtis)

    total_interval = Interval{T, Closed, Open}[Interval{T, Closed, Open}(start_time, stop_time)]
    total_interval_set = IntervalSet(total_interval)

    gti_interval = Interval{T, Closed, Open}[]
    for ig in eachrow(gtis)
        push!(gti_interval,Interval{T, Closed, Open}(ig[1],ig[2]))
    end
    gti_interval_set = IntervalSet(gti_interval)

    bti_interval_set = setdiff(total_interval_set, gti_interval_set)

    btis = Vector{T}[]

    for interval in bti_interval_set.items
        push!(btis, [first(interval), last(interval)])
    end

    return mapreduce(permutedims, vcat, btis)
end

function time_intervals_from_gtis(gtis::AbstractMatrix{<:Real}, segment_size::Real;
                                  fraction_step::Real=1, epsilon::Real=1e-5)  
    spectrum_start_times = Float64[]

    gti_low = @view gtis[:,1]
    gti_up = @view gtis[:,2]

    for (g1,g2) in zip(gti_low,gti_up)
        if g2 - g1 + epsilon < segment_size
            continue
        end

        newtimes = range(g1, g2 - segment_size + epsilon, step = segment_size* fraction_step)
        append!(spectrum_start_times,newtimes)
    end
    return spectrum_start_times, spectrum_start_times .+ segment_size
end

function calculate_segment_bin_start(startbin::Integer, stopbin::Integer,
                                     nbin::Integer; fraction_step::Real=1)
    st = floor.(range(startbin, stopbin, step=Int(nbin * fraction_step)))
    if st[end] == stopbin
        pop!(st)
    end
    if st[end] + nbin > stopbin
        pop!(st)
    end
    return st
end

function bin_intervals_from_gtis(gtis::AbstractMatrix{<:Real}, segment_size::Real,
                                 time::AbstractVector{<:Real}; dt=nothing, 
                                 fraction_step::Real=1, epsilon::Real=0.001)
    if isnothing(dt)
        dt = Statistics.median(diff(time))
    end

    epsilon_times_dt = epsilon * dt
    nbin = round(Int, segment_size / dt)

    spectrum_start_bins = Int[]

    gti_low = @view(gtis[:, 1]) .+ (dt ./ 2 .- epsilon_times_dt)
    gti_up = @view(gtis[:, 2]) .- (dt ./ 2 .- epsilon_times_dt)

    for (g0, g1) in zip(gti_low, gti_up)
        if (g1 - g0 .+ (dt + epsilon_times_dt)) < segment_size
            continue
        end
        startbin, stopbin = searchsortedfirst.(Ref(time), [g0, g1])
        startbin -= 1
        if stopbin > length(time)
            stopbin = length(time)
        end

        if time[startbin+1] < g0
            startbin += 1
        end
        # Would be g[1] - dt/2, but stopbin is the end of an interval
        # so one has to add one bin
        if time[stopbin] > g1
            stopbin -= 1
        end

        newbins = calculate_segment_bin_start(
            startbin, stopbin, nbin, fraction_step=fraction_step)
        
        append!(spectrum_start_bins,newbins)
    end 
    return spectrum_start_bins, spectrum_start_bins.+nbin 
end

@resumable function generate_indices_of_segment_boundaries_unbinned(times::AbstractVector{<:Real},
                                                                    gti::AbstractMatrix{<:Real},
                                                                    segment_size::Real)
    start, stop = time_intervals_from_gtis(gti, segment_size)

    startidx = searchsortedfirst.(Ref(times), start)
    stopidx = searchsortedfirst.(Ref(times), stop)

    for (s, e, idx0, idx1) in zip(start, stop, startidx, stopidx)
        @yield s, e, idx0, idx1
    end
end

@resumable function generate_indices_of_segment_boundaries_binned(times::AbstractVector{<:Real},
                                                                  gti::AbstractMatrix{<:Real},
                                                                  segment_size::Real; dt=nothing)
    startidx, stopidx = bin_intervals_from_gtis(gti, segment_size, times;
                                                dt=dt)

    if isnothing(dt)
        dt = 0
    end
    for (idx0, idx1) in zip(startidx, stopidx)
        @yield times[idx0+1] - dt / 2, times[min(idx1, length(times) - 1)] - dt / 2,idx0, idx1
    end
end

#gti interface

"""
    apply_gtis(el::EventList{T}, gtis::AbstractMatrix{<:Real}) where T

Apply Good Time Intervals (GTIs) to an EventList, returning an array of EventLists - one for each GTI.
"""
function apply_gtis(el::EventList{T}, gtis::AbstractMatrix{<:Real}) where T
    check_gtis(gtis)
    
    result = EventList{T}[]
    
    for i in 1:size(gtis, 1)
        gti_start, gti_end = gtis[i, 1], gtis[i, 2]
        mask = (el.times .>= gti_start) .& (el.times .<= gti_end)
        
        filtered_times = el.times[mask]
        filtered_energies = el.energies[mask]
        
        push!(result, EventList{T}(
            el.filename,
            filtered_times,
            filtered_energies,
            deepcopy(el.metadata))
        )
    end
    
    return result
end

"""
    apply_gtis(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real}) where T

Apply Good Time Intervals (GTIs) to a LightCurve, returning an array of LightCurves - one for each GTI.
"""
function apply_gtis(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real}) where T
    check_gtis(gtis)
    
    result = LightCurve{T}[]
    
    for i in 1:size(gtis, 1)
        gti_start, gti_end = gtis[i, 1], gtis[i, 2]
        mask = (lc.timebins .>= gti_start) .& (lc.timebins .<= gti_end)
        
        filtered_timebins = lc.timebins[mask]
        filtered_counts = lc.counts[mask]
        filtered_errors = lc.count_error[mask]
        
        push!(result, LightCurve{T}(
            filtered_timebins,
            filtered_counts,
            filtered_errors,
            lc.err_method)
        )
    end
    
    return result
end

"""
    fill_bad_time_intervals!(el::EventList{T}, gtis::AbstractMatrix{<:Real};
                           fill_value::T=zero(T)) where T

Fill bad time intervals (outside GTIs) in an EventList with a fill value.
"""
function fill_bad_time_intervals!(el::EventList{T}, gtis::AbstractMatrix{<:Real};
                                 fill_value::T=zero(T)) where T
    check_gtis(gtis)
    
    outside_mask = fill(true, length(el.times))
    for i in 1:size(gtis, 1)
        gti_start, gti_end = gtis[i, 1], gtis[i, 2]
        outside_mask .&= .!((el.times .>= gti_start) .& (el.times .<= gti_end))
    end
    
    el.energies[outside_mask] .= fill_value
    
    return nothing
end

"""
    fill_bad_time_intervals!(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real};
                           fill_value::T=zero(T)) where T

Fill bad time intervals (outside GTIs) in a LightCurve with a fill value.
"""
function fill_bad_time_intervals!(lc::LightCurve{T}, gtis::AbstractMatrix{<:Real};
                                 fill_value::T=zero(T)) where T
    check_gtis(gtis)
    
    outside_mask = fill(true, length(lc.timebins))
    for i in 1:size(gtis, 1)
        gti_start, gti_end = gtis[i, 1], gtis[i, 2]
        outside_mask .&= .!((lc.timebins .>= gti_start) .& (lc.timebins .<= gti_end))
    end
    
    lc.counts[outside_mask] .= fill_value
    lc.count_error[outside_mask] .= zero(T)
    
    return nothing
end