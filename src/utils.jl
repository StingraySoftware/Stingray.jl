import Base

function sum_if_not_none_or_initialize(A,B)
    if isnothing(A) 
        return deepcopy((B))
    end
    return A + B
end

sqrt(x::Real) = x < 0.0 ? NaN : Base.sqrt(x)
