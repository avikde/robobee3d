
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
εs = [1, 0.005, 0.001] # IC, dyn, symm
state0 = traj0, param0
state1 = cu.optboth(nothing, m, opt, state0..., εs; step=1e2)
state2 = cu.optboth(nothing, m, opt, state1..., εs; step=1e2)
state3 = cu.optboth(nothing, m, opt, state2..., εs; step=1e2)

states = [state0, state1, state2, state3]
# Results ---

trajs = first.(states)
params = last.(states)

println("Objectives: ", [(_ro = cu.robj(m, opt, trajs[i], params[i]); _ro'*_ro) for i = 1:length(trajs)])
println("Params: ", params)

pl1 = plotTrajs(m, opt, trajt, params, trajs)
pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

l = @layout [grid(2,1) a]
plot(pl1..., pl2, layout=l, size=(900,400))
