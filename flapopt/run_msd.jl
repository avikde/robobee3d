
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel(15, 0.1, 1, 50)
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, 0.2, 1, :symmetric, 1e-8, false)
# opt = cu.OptOptions(false, 0.2, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == :symmetric ? 17 : 34
param0 = [1.0] # k

trajt, traj0, trajt1 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0)

# Add a priority matrix for weighting the constraints in the param opt step
nc = ny * (N+2)
Q = zeros(1,nc)
for i = 1:N
	Q[1, i*ny + 1:(i+1)*ny] = [0, 1]
end

# IPOPT
εs = [1, 0.005, 0.001] # IC, dyn, symm
# state0 = traj0, param0
# state1 = cu.optboth(nothing, m, opt, state0..., εs, Q; step=1e2)
# states = [state0, state1]
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs, Q; step=1e1))
# m.umax = 25.0
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs, Q; step=1e1))
# # m.umax = 20.0
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs, Q; step=1e1))
# push!(states, cu.optboth(nothing, m, opt, states[end]..., εs, Q; step=1e1))

# # Hand-craft a trajectory+param
# paramHC = [resFreq(m, opt, traj0)]
# prob = cu.ipoptsolve(m, opt, first(states[end]), paramHC, εs, :traj; print_level=1, nlp_scaling_method="none")
# push!(states, (prob.x, paramHC))

# # # Numerically try: fix param, test for IC, set controller=0 and look @traj
# # trajOL = createOLTraj(m, opt, traj0, paramHC)
# # push!(states, (trajOL, paramHC))
# # Results ---

# trajs = first.(states)
# paramss = last.(states)

# println("Objectives: ", [(_ro = cu.robj(m, opt, trajs[i], paramss[i]); _ro'*_ro) for i = 1:length(trajs)])
# println("Params: ", paramss)

# pl1 = plotTrajs(m, opt, trajt, paramss, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, paramss, trajs)

# l = @layout [grid(2,1) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# Traj
Hh = hcat([Hdes(m, 0.1, t) for t in trajt1]...)
plot(trajt, Hh', marker=:auto)

# # Test naive
# ktest, os1 = cu.optnaive(nothing, m, opt, traj0, εs)
# kr1 = resStiff(m, opt, traj0)
# m.bσ = 8
# ktest, os3 = cu.optnaive(nothing, m, opt, traj0, εs)
# pl1 = plot(ktest, [os1, os3], marker=:auto, xlabel="p", ylabel="obj")
# vline!(pl1, [kr1])
# plot(pl1)
