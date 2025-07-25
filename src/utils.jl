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
    merge_overlapping_gtis(gtis::Matrix{Float64}) -> Matrix{Float64}

Merge overlapping or touching GTI intervals into continuous segments.

Takes a matrix of GTI intervals and combines any that overlap or touch into single
continuous intervals. The input GTIs are automatically sorted by start time.

# Arguments
- `gtis::Matrix{Float64}`: Matrix of GTI boundaries where each row is [start_time, stop_time]

# Returns
- `Matrix{Float64}`: Matrix with merged intervals, sorted by start time

# Examples
```julia
# Basic merging
gtis = [0.0 2.0; 1.0 3.0; 4.0 5.0]
result = merge_overlapping_gtis(gtis)
# Returns: [0.0 3.0; 4.0 5.0]

# Touching intervals
gtis = [0.0 1.0; 1.0 2.0; 2.0 3.0]
result = merge_overlapping_gtis(gtis)
# Returns: [0.0 3.0]
```
"""
function merge_overlapping_gtis(gtis::Matrix{Float64})
    if size(gtis, 1) <= 1
        return gtis
    end
    
    # Sort by start time
    sort_indices = sortperm(view(gtis, :, 1))
    sorted_gtis = gtis[sort_indices, :]
    
    # Pre-allocate merged array
    merged = Matrix{Float64}(undef, size(gtis, 1), 2)
    merged_count = 0
    
    current_start = sorted_gtis[1, 1]
    current_stop = sorted_gtis[1, 2]
    
    for i in 2:size(sorted_gtis, 1)
        start_time = sorted_gtis[i, 1]
        stop_time = sorted_gtis[i, 2]
        
        # Check if intervals overlap or touch (with small tolerance)
        if start_time <= current_stop + 1e-6
            current_stop = max(current_stop, stop_time)
        else
            # No overlap, add current interval to merged
            merged_count += 1
            merged[merged_count, 1] = current_start
            merged[merged_count, 2] = current_stop
            current_start = start_time
            current_stop = stop_time
        end
    end
    
    # Add the final interval
    merged_count += 1
    merged[merged_count, 1] = current_start
    merged[merged_count, 2] = current_stop
    
    return merged[1:merged_count, :]
end
