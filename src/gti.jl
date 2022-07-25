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
