"""
    qvalues(P, λ = 0:01:0.01:0.99)

Calculate q-values from a set of p-values P. λ can be a range of values to test,
or a fixed value.
"""
function qvalues(P, λ = 0:01:0.01:0.99)
    order = sortperm(P)
    P = P[order]
    m = length(P)
    pm = P[m]
    # π̂ = sum(P .> λ) / (m * (1 - λ))
    π̂₀ = estimate_π̂₀(P, λ)
    q̂_pm = π̂ * pm
    Q̂ = [q̂_pm]
    for i in (m-1):-1:1
        q̂ = min(
            π̂ * m * P[i] / i,
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
π̂₀(P, λ) = mean(P .>= λ) / (1 - λ)


"""
    estimate_π̂₀(P, λs = 0.01:0.01:0.99)

Estimate the proportion of true null hypotheses, π̂₀, given a range of λs.
"""
function estimate_π̂₀(P, λs = 0.01:0.01:0.99)
    π̂s = Float64[]
    for λ in λs
        # push!(π̂s, min(1, π̂₀(P, λ)))
        push!(π̂s, π̂₀(P, λ))
    end
    # return π̂s
    itp = interpolate(π̂s, BSpline(Cubic(Natural())), OnGrid())
    sitp = scale(itp, λs)
    return sitp
end
