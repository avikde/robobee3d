
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
opt = cu.OptOptions(false, 0.2, 1, :symmetric, 1e-8, false)
# opt = cu.OptOptions(false, 0.2, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == :symmetric ? 17 : 34
param0 = [2.0, 20.0] # cbar, T

trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], params0)

# # 
# trajei = eulerIntegrate(m, opt, traj0, params0)
# pl1 = plotTrajs(m, opt, trajt, params0, traj0, trajei)

# traj opt ------------------------------------

# IPOPT
εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# traj1 = cu.ipoptsolve(m, opt, traj0, params0, εs, :traj)
# trajs = [traj0, traj1]
# params = [param0, param0]

# naive param opt
εsp = [100.0, 0.005, 100.0] # IC, dyn, symm
param1 = cu.ipoptsolve(m, opt, traj0, param0, εs, :param)
trajs = [traj0, traj0]
params = [param0, param1]

# # with Coros g-preferred param opt
# param1 = cu.paramopt(m, opt, traj1, param0, εs; step=0.3)
# traj2 = cu.ipoptsolve(m, opt, traj1, param1, εs, :traj)

# trajs = [traj0, traj1, traj2]
# params = [param0, param0, param1]

# # Custom solver
# wkt = cu.OptWorkspace(cu.Ntraj(m, opt, N), (N+2)*ny)
# @time traj2 = cu.csSolve!(wkt, m, opt, traj0, params0, :traj; Ninner=30, μs=[1e5])
# trajs = [traj0, traj1, traj2]
# params = [params0, params0, params0]
# # trajs, params, wkt = cu.csAlternateSolve(m, opt, traj0, params0, 1; μst=[1e6], Ninnert=30, μsp=[1e-2,1e-2], Ninnerp=2)

# pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# display(params)

# paramopt -------------------------------------


# visualize --------------------------------------------

println("Objectives: ", [(_ro = cu.robj(m, opt, tt, params0); _ro ⋅ _ro) for tt in trajs])
println("Params: ", params)

# animateTrajs(m, opt, params, trajs)
pl1 = plotTrajs(m, opt, trajt, params, trajs)
pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

l = @layout [grid(2,2) a]
plot(pl1..., pl2, layout=l, size=(900,400))

# plot(pl1...)
# gui()
