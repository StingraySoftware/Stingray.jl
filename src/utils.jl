function sum_if_not_none_or_initialize(A,B)
    if isnothing(A) 
        return deepcopy(copy((B)))
    end
    return A + B
end
