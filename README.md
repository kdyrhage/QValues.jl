# QValues.jl
## Installation
QValues.jl can be install from within Julia using
```julia
Pkg.clone("git@github.com:kdyrhage/QValues.jl.git")
```

## Usage
Examples usage:
```julia
using QValues
using RDatasets

# Load a dataset that contains p-values. This set actually
# contains values 0 < p < 0.05, so multiply by 20 to get
# 0 < p < 1
manhattan = dataset("gap", "mhtdata")
P = 20 * manhattan[:P]

Q = qvalues(P, 0.05:0.05:0.95, method = :bootstrap)

using Plots
scatter(P, Q, xlabel = "p-value", ylabel = "q-value", marker = (2, :black))
```

Three methods are available for estimating the null hypothesis rate: ```:spline```, ```:bootstrap```, and ```:fast```.
