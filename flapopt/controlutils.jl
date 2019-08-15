
module controlutils

using LinearAlgebra
using ForwardDiff
import Ipopt # keep in its namespace
# using PositiveFactorizations
using Plots

include("Model.jl")
include("OptIPOPT.jl")
include("OptCustom.jl")
include("OptParams.jl")

end
