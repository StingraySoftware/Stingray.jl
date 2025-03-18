# this function is returning a scalar sum instead of maintaining the array structure.
# this is to prevent array interface problem
function sum_if_not_none_or_initialize(current, new_value)
    if isnothing(new_value)
        return current  # If new_value is nothing, just return current
    end

    if isnothing(current)
        return copy(new_value)  # Initialize if current is nothing
    end

    if isa(current, AbstractArray) && isa(new_value, AbstractArray)
        return current .+ new_value  # Element-wise addition for arrays
    elseif isa(current, Number) && isa(new_value, Number)
        return current + new_value  # Simple number addition
    else
        error("sum_if_not_none_or_initialize: Type mismatch between current=$(typeof(current)) and new_value=$(typeof(new_value))")
    end
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
        push!(idx, condition.size + 1)
    end

    # Reshape the result into two columns
    return reshape(idx, 2, length(idx) ÷ 2)'
end
