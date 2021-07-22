# Averages the xs and associated ys in each grouping of x
function bin_data(x, y; Δ=0.1)
    idxs = idxs_of_groups(x, Δ)
    groups = [(x[idx], y[idx]) for idx in idxs]
    x_avg, y_avg = Float64[], Measurement{Float64}[]
    for group in groups
        push!(x_avg, mean(group[1]))
        push!(y_avg, weightedmean(group[2]))
    end

    return x_avg, y_avg
end

# Given a sorted vector `v`, return the set of groups (by idx) satisfying the condition that all points in that group are ≤ Δv of the first point
function idxs_of_groups(v, Δv)
    idx_start = 1
    groups = Vector{Int}[]
    group = [idx_start]
    for i in 2:length(v)
        if v[i] - v[idx_start] ≤ Δv
            push!(group, i)
        else
            push!(groups, group)
            idx_start = i
            group = [i]
        end
    end
    push!(groups, group)

    return groups
end
