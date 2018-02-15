"""
    qvalues(P, λ = 0.05:0.01:0.95)

Calculate q-values from a set of p-values P. λ can be a range of values to test,
or a fixed value.
"""
function qvalues(P, λ = 0.05:0.01:0.99, π̂₀ = 0.0;
                 method = :bootstrap, B = 100, fraction = 0.6, γ = 0.05)
    order = sortperm(P)
    P = P[order]
    m = length(P)
    if π̂₀ == 0.0
        if length(λ) == 1
            π̂₀ = π̂(P, λ)
        else
            if method == :bootstrap
                π̂₀ = bootstrap_π̂₀(P, γ, λ; B = B, f = fraction)
            elseif method == :fast
                π̂₀ = fast_π̂₀(P, λ)
            elseif method == :spline
                π̂₀ = spline_π̂₀(P, λ)
            else
                error("'$(String(method))' is not a valid method for estimating π̂₀")
            end
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
    spline_π̂₀(P, λs = 0.05:0.01:0.95)

Estimate the proportion of true null hypotheses, π̂₀, given a range of λs.
"""
function spline_π̂₀(P, λs = 0.05:0.01:0.95)
    π̂s = Float64[]
    for λ in λs
        push!(π̂s, π̂(P, λ))
    end
    spl = fit(SmoothingSpline, collect(λs), π̂s, 0.015) # 0.015 might not be a good λ for every dataset
    return SmoothingSplines.predict(spl, 1.0)
end


"""
    pFDRλ(P, λ, γ)

Estimate pFDR(γ) for a given λ.
"""
function pFDRλ(P, λ, γ)
    m = length(P)
    π̂₀ = π̂(P, λ)
    P̂rPγ = max(sum(P .<= γ), 1.) / m
    pF̂DRλ = (π̂₀ * γ) / (P̂rPγ * (1 - (1 - γ)^m))
    return min(1., pF̂DRλ)
end


"""
    bootstrap_π̂₀(P, γ = 0.05, λs = 0.05:0.01:0.95; B = 100, f = 0.6)

Use bootstrap method to estimating the best λ.

Keyword arguments `B` controls the number of bootstraps to run for each λ, while
`f` is the fraction of the original dataset to use for each bootstrap.
"""
function bootstrap_π̂₀(P, γ = 0.05, λs = 0.01:0.01:0.95; B = 100, f = 0.6)
    nsamples = floor(Int, f * length(P))
    Pb = zeros(Float64, nsamples)
    pF̂DRλ = Float64[]
    for λ in λs
        push!(pF̂DRλ, pFDRλ(P, λ, γ))
    end
    M̂SEλ = Float64[]
    for λ in λs
        pF̂DRλb = Float64[]
        for b in 1:B
            StatsBase.seqsample_a!(P, Pb)
            push!(pF̂DRλb, pFDRλ(Pb, λ, γ))
        end
        push!(M̂SEλ, (1/B) * sum((pF̂DRλb .- minimum(pF̂DRλ)).^2))
    end
    λ̂ = minimum(λs[indmin(M̂SEλ)])
    return π̂(P, λ̂)
end


"""
    fast_π̂₀(P, λs)

Estimate π̂₀ using a fast method.

This is the "bootstrap" method implemented in Storey's qvalue package in R. It
doesn't actually perform bootstrapping, but it is very fast.
"""
function fast_π̂₀(P, λs = 0.01:0.01:0.95)
    m = length(P)
    π̂s = [π̂(P, λ) for λ in λs]
    minπ̂ = quantile(π̂s, 0.1)
    W = [sum(P .>= λ) for λ in λs]
    mse = (W ./ (m^2 .* (1 .- λs).^2)) .* (1 .- W ./ m) .+ (π̂s .- minπ̂).^2
    π̂₀ = min(π̂s[indmin(mse)]..., 1)
end
