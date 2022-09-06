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

function rebin_data(x::AbstractVector{<:Real}, y::AbstractVector{<:Integer}, dx_new::Real;
                    y_err::AbstractVector{<:Real}=Float64[], method::String = "sum", dx::Real=0)
    if isempty(y_err)
        y_err = zero(y)
    end

    if dx isa AbstractVector
        dx_old = dx
    elseif dx!=0
        dx_old = [dx]
    else
        dx_old = diff(x)
    end

    if any(dx_new .< dx_old)
        throw(ArgumentError("New frequency resolution must be larger than old frequency resolution."))
    end

    # left and right bin edges
    # assumes that the points given in `x` correspond to
    # the left bin edges
    xedges = vcat(x, x[end]+dx_old[end])

    # new regularly binned resolution
    xbin = range(xedges[begin], xedges[end], step = dx_new)

    output = zeros(length(xbin) - 1)
    outputerr = zeros(length(xbin) - 1)
    step_size = zeros(length(xbin) - 1)

    all_x = searchsortedfirst.(Ref(xedges), xbin)
    min_inds = @view(all_x[begin:end-1])
    max_inds = @view(all_x[begin+1:end])
    xmins = @view(xbin[begin:end-1])
    xmaxs = @view(xbin[begin+1:end])

    for (i, (xmin, xmax, min_ind, max_ind)) in enumerate(zip(xmins, xmaxs, min_inds, max_inds))
        filtered_y = @view(y[min_ind:max_ind-1])
        filtered_yerr = @view(y_err[min_ind:max_ind-1])
        output[i] = sum(filtered_y)
        outputerr[i] = sum(filtered_yerr)
        step_size[i] = max_ind - 1 - min_ind

        prev_dx = xedges[min_ind+1] - xedges[min_ind]
        prev_frac = (xedges[min_ind] - xmin)/prev_dx
        output[i] += y[min_ind]*prev_frac
        outputerr[i] += y_err[min_ind]*prev_frac
        step_size[i] += prev_frac

        if !(max_ind == length(xedges))
            dx_post = xedges[max_ind] - xedges[max_ind-1]
            post_frac = (xmax-xedges[max_ind-1])/dx_post
            output[i] += y[max_ind-1]*post_frac
            outputerr[i] += y_err[max_ind-1]*post_frac
            step_size[i] += post_frac
        end
    end

    if method in ["mean", "avg", "average", "arithmetic mean"]
        ybin = output / step_size
        ybinerr = sqrt.(outputerr) / step_size
    elseif method == "sum"
        ybin = output
        ybinerr = sqrt.(outputerr)
    else
        throw(ArgumentError("Method for summing or averaging not recognized. Please enter either 'sum' or 'mean'."))
    end
    tseg = x[end] - x[begin] + dx_old[end]

    if (tseg / dx_new % 1) > 0
        pop!(ybin)
        pop!(ybinerr)
        pop!(step_size)
    end

    dx_var = Statistics.var(dx_old) / Statistics.mean(dx_old)

    if size(dx_old) == 1 || dx_var < 1e-6
        new_step = step_size[1]
    end

    new_x0 = (x[1] - (0.5 * dx_old[1])) + (0.5 * dx_new)
    xbin = 1:length(ybin) * dx_new .+ new_x0

    return xbin, ybin, ybinerr, new_step
end
