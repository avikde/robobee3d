
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
using SparseArrays, OSQP, LinearAlgebra # temp
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

fdes = 0.15
trajt, traj0, trajt1 = createInitialTraj(m, opt, N, fdes, [1e3, 1e2], param0)

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

# Optimize params directly for traj
Hh = [Hdes(m, fdes, t) for t in trajt1]
np = 3
P1 = zeros(np,np)
for Hk in Hh
	P1 .= P1 + Hk * Hk'
end
mo = OSQP.Model()
OSQP.setup!(mo; P=sparse(P1), q=zeros(np), A=sparse(1:np, 1:np, ones(np)), l=[0; 0; m.mb], u=[Inf; m.bσ; m.mb])
res = OSQP.solve!(mo)
resf = sqrt(res.x[1]/m.mb)/(2*π)
println("Res freq = ", resf)

# Plot the quadratic
size = 100
x = range(0, stop=20, length=size)
y = range(0, stop=10, length=size)
ff(x,y) = [x;y;m.mb]' * P1 * [x;y;m.mb]
pl1 = contour(x, y, ff)
vline!(pl1, [res.x[1]])

# plot(trajt, Hh', marker=:auto)

# # Test naive
# ktest, os1 = cu.optnaive(nothing, m, opt, traj0, εs)
# kr1 = resStiff(m, opt, traj0)
# m.bσ = 8
# ktest, os3 = cu.optnaive(nothing, m, opt, traj0, εs)
# pl1 = plot(ktest, [os1, os3], marker=:auto, xlabel="p", ylabel="obj")
# vline!(pl1, [kr1])
# plot(pl1)
