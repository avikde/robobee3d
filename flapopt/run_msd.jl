
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel()
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, 0.2, 1, :symmetric, 1e-8, false)
# opt = cu.OptOptions(false, 0.2, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == :symmetric ? 17 : 34
param0 = [1.0] # k

trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0)


# IPOPT
εs = [0.05, 0.005, 0.001] # IC, dyn, symm
prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1, nlp_scaling_method="none")
traj1 = prob.x

trajs = [traj0, traj1]
params = [param0, param0]
pl1 = plotTrajs(m, opt, trajt, params, trajs)

