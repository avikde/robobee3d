
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("Wing2DOF.jl")

# create an instance
m = Wing2DOFModel()
ny, nu = cu.dims(m)
opt = cu.OptOptions(true, 0.1, 3, :symmetric, 1e-8, false)
# opt = cu.OptOptions(false, 0.2, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == :symmetric ? 11 : 22
params0 = [2.0, 20.0] # cbar, T

trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], params0)

# 
trajei = eulerIntegrate(m, opt, traj0, params0)
pl1 = plotTrajs(m, opt, trajt, params0, traj0, trajei)

# # IPOPT
# prob = cu.nloptsetup(m, opt, traj0, params0)
# status = cu.nloptsolve(prob)
# trajs = [traj0, prob.x]
# animateTrajs(m, opt, params0, traj0, prob.x)
# pl1 = plotTrajs(m, opt, trajt, params0, traj0, prob.x)
# pl2 = cu.visualizeConstraintViolations(m, opt, params0, traj0, prob.x)

# # Custom solver
# trajs, params, wkt = cu.csAlternateSolve(m, opt, traj0, params0, 1; μst=[1e6, 1e3], Ninnert=10, μsp=[1e-2,1e-2], Ninnerp=2)

# animateTrajs(m, opt, params0, [view(trajs, :, i) for i in 1:size(trajs, 2)]...)

# pl1 = plotTrajs(m, opt, trajt, params0, (trajs[:,i] for i = 1:size(trajs,2))...)
# # pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# # display(params)

# visualize constraint violations

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

plot(pl1...)
# gui()
