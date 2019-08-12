
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using Plots, BenchmarkTools
import Ipopt
gr() # backend
using Revise # while developing
import controlutils
cu = controlutils
includet("Wing2DOF.jl")

# create an instance
m = Wing2DOFModel()
ny, nu = cu.dims(m)
opt = cu.OptOptions(true, 0.1, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == cu.SYMMETRIC ? 11 : 22
params0 = [2.0, 20.0] # cbar, T

trajt, traj0 = createInitialTraj(m, N, 0.15, [1e3, 1e2], params0)

# wkt = cu.OptWorkspace((N+1)*ny + N*nu + 1, (N+2)*ny)
# cu.csSolve!(wk, m, opt, traj0, params0, cu.WRT_TRAJ)

# IPOPT
prob = cu.nloptsetup(m, opt, traj0, params0; tol=1e-2, constr_viol_tol=1e1, acceptable_constr_viol_tol=1e1, acceptable_dual_inf_tol=1e1)
status = cu.nloptsolve(prob)
trajs = [traj0, prob.x]
pl1 = plotTrajs(m, opt, trajt, params0, traj0, prob.x)
pl2 = cu.visualizeConstraintViolations(m, opt, params0, traj0, prob.x)

# # Custom solver
# trajs, params, wkt = cu.csAlternateSolve(m, opt, traj0, params0, 1; μst=[1e6, 1e3], Ninnert=10, μsp=[1e-2,1e-2], Ninnerp=2)

# pl1 = plotTrajs(m, opt, trajt, params0, (trajs[:,i] for i = 1:size(trajs,2))...)
# # pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# # display(params)

# visualize constraint violations

l = @layout [grid(2,2) a]
plot(pl1..., pl2, layout=l, size=(900,400))

# plot(pl1...)
gui()
