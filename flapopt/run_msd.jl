
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
states = [state0, state1]
push!(states, cu.optboth(nothing, m, opt, states[end]..., εs; step=1e1))
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs; step=1e1))
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs; step=1e1))
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs; step=1e1))

# Hand-craft a trajectory+param
paramHC = [4*resFreq(m, opt, traj0)]
prob = cu.ipoptsolve(m, opt, first(states[end]), paramHC, εs, :traj; print_level=1, nlp_scaling_method="none")
# push!(states, (prob.x, paramHC))

# Numerically try: fix param, test for IC, set controller=0 and look @traj
trajOL = createOLTraj(m, opt, traj0, paramHC)
push!(states, (trajOL, paramHC))
# Results ---

trajs = first.(states)
paramss = last.(states)

println("Objectives: ", [(_ro = cu.robj(m, opt, trajs[i], paramss[i]); _ro'*_ro) for i = 1:length(trajs)])
println("Params: ", paramss)

pl1 = plotTrajs(m, opt, trajt, paramss, trajs)
pl2 = cu.visualizeConstraintViolations(m, opt, paramss, trajs)

l = @layout [grid(2,1) a]
plot(pl1..., pl2, layout=l, size=(900,400))
