function sum_if_not_none_or_initialize(A,B)
    if isnothing(A) 
        return deepcopy((B))
    end
    return A + B
end

function contiguous_regions(condition::AbstractVector{Bool})
    # Find the indicies of changes in "condition"
    d = diff(condition)
    idx = findall(!iszero, d)

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx .+= 1

    if condition[1]
        # If the start of condition is True prepend a 0
        pushfirst!(idx, 1)
    end

    if condition[end]
        # If the end of condition is True, append the length of the array
        push!(idx, condition.size+1)
    end

    # Reshape the result into two columns
    return reshape(idx,2,length(idx) ÷ 2)'
end
#UTIL function used in recpies
"""
    merge_overlapping_gtis(gtis::Matrix{Float64})::Matrix{Float64}
   
Merge overlapping Good Time Intervals (GTIs).
This function takes a matrix of GTIs, where each row represents a start and stop time,
and returns a new matrix with non-overlapping GTIs.
The input matrix `gtis` should have two columns, where the first column is the start time and the second column is the stop time.
If the input matrix has one or zero rows, it is returned unchanged.
The function sorts the GTIs by start time, then iterates through them, merging any overlapping intervals.
The output is a matrix with the same two-column format, containing the merged GTIs.
"""
function merge_overlapping_gtis(gtis::Matrix{Float64})::Matrix{Float64}
    if size(gtis, 1) <= 1
        return gtis
    end
    sort_indices = sortperm(gtis[:, 1])
    sorted_gtis = gtis[sort_indices, :]
    
    merged = Matrix{Float64}(undef, 0, 2)
    current_start = sorted_gtis[1, 1]
    current_stop = sorted_gtis[1, 2]
    
    for i in 2:size(sorted_gtis, 1)
        start_time = sorted_gtis[i, 1]
        stop_time = sorted_gtis[i, 2]
        if start_time <= current_stop + 1e-6  # Small tolerance for floating point
            # Merge intervals - extend the current stop time
            current_stop = max(current_stop, stop_time)
        else
            merged = vcat(merged, [current_start current_stop])
            current_start = start_time
            current_stop = stop_time
        end
    end
    
    merged = vcat(merged, [current_start current_stop])
    
    return merged
end