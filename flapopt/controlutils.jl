
module controlutils

using LinearAlgebra
using ForwardDiff
import Ipopt # keep in its namespace
# using PositiveFactorizations
using Plots; gr()

include("Model.jl")
include("ModelOptimizer.jl")

end
