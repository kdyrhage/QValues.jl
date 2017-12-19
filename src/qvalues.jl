"""
    qvalues(P, λ = 0.05:0.01:0.95)

Calculate q-values from a set of p-values P. λ can be a range of values to test,
or a fixed value.
"""
function qvalues(P, λ = 0.05:0.01:0.95, π̂₀ = 0.0)
    order = sortperm(P)
    P = P[order]
    m = length(P)
    if π̂₀ == 0.0
        if length(λ) == 1
            π̂₀ = π̂(P, λ)
        else
            π̂₀ = estimate_π̂₀(P, λ)
        end
    end
    q̂ₚₘ = π̂₀ * P[m]
    Q̂ = [q̂ₚₘ]
    for i in (m-1):-1:1
        q̂ = min(
            π̂₀ * m * P[i] / i,
            Q̂[end]
        )
        push!(Q̂, q̂)
    end
    return Q̂[sortperm(reverse(order))]
end


"""
    π̂₀(P, λ)

Calculate π̂₀ for a vector of p-values P at a given λ.
"""
π̂(P, λ) = mean(P .> λ) / (1 - λ)


"""
    estimate_π̂₀(P, λs = 0.05:0.01:0.95)

Estimate the proportion of true null hypotheses, π̂₀, given a range of λs.
"""
function estimate_π̂₀(P, λs = 0.05:0.01:0.95)
    π̂s = Float64[]
    for λ in λs
        push!(π̂s, π̂(P, λ))
    end
    spl = fit(SmoothingSpline, collect(λs), π̂s, 0.015) # 0.015 might not be a good λ for every dataset
    return predict(spl, 1.0)
end
