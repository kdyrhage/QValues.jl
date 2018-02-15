# QValues.jl
## Installation
QValues and its one dependency can be installed from within Julia using
```julia
Pkg.add("SmoothingSplines")
Pkg.clone("git@github.com:kdyrhage/QValues.jl.git")
```

## Usage
QValues only exports one function:
```julia
# Calculate q-values given a set of p-values P:
qvalues(P)

# A range or vector containing lambda values to be tested can also be included:
qvalues(P, 0.01:0.01:0.99)

# The method for estimating the rate of true nulls can be changed, and for
# :bootstrap additional variables can be set:
qvalues(P, method = :bootstrap, B = 1000, fraction = 0.8)

# Other estimation methods ignore these options, so the following two are equivalent:
qvalues(P, method = :fast, B = 1, fraction = 0.5)
qvalues(P, method = :fast)
```

Three methods are available for estimating the null hypothesis rate: ```:spline```, ```:bootstrap```, and ```:fast```. ```:spline``` and ```:bootstrap``` use the methods described in [1] and [2], respectively, while ```:fast``` uses the "bootstrap" method used in the [qvalue package for R](https://github.com/StoreyLab/qvalue/) (which does not actually involve any bootstrapping). ```:bootstrap``` is much slower than the other two methods. The number of bootstraps can be adjusted with the ```B``` keyword, and the size of the fraction used for each iteration can be adjusted with ```fraction```.

Examples usage:
```julia
using QValues
using RDatasets
using Plots

# Load a dataset that contains p-values. This set actually contains
# values in 0 < p < 0.05, so multiply by 20 to get 0 < p < 1
manhattan = dataset("gap", "mhtdata")
P = 20 * manhattan[:P]

Q = qvalues(P, 0.05:0.05:0.95, method = :bootstrap)

scatter(P, Q, xlabel = "p-value", ylabel = "q-value", marker = (2, :black))
```

Example usage without external packages:
```julia
using QValues

# Rough simulation of p-values for true null (p0) and true positives (p1)
P0 = rand(1000)
P1 = vcat([rand(div(1000, 2i)) / i^2 for i in 4:9]...)
P = [P0; P1]

Q = qvalues(P)
sum(Q .<= 0.05)
```
