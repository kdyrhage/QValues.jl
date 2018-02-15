using QValues
using Base.Test

# true null
p0 = collect(0.01:0.01:1.0)
# true positive
p1 = collect(0.001:0.001:0.05)
p = vcat(p0, p1)

# rate of true null
π₀ = length(p0) / length(p)

@test QValues.spline_π̂₀(p) ≈ π₀
@test QValues.bootstrap_π̂₀(p) ≈ π₀
@test QValues.fast_π̂₀(p) ≈ π₀
