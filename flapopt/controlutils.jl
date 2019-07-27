
module controlutils

using LinearAlgebra
using ForwardDiff
import Ipopt # keep in its namespace

include("Model.jl")
include("ModelOptimizer.jl")

end
