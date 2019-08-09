
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
opt = cu.OptOptions(true, 0.1, 1, cu.SYMMETRIC, 1e-4)
N = opt.boundaryConstraint == cu.SYMMETRIC ? 11 : 22
params0 = [2.0, 20.0] # cbar, T

trajt, traj0 = createInitialTraj(m, N, 0.15, [1e3, 1e2], params0)

wk = cu.OptWorkspace((N+1)*ny + N*nu + 1, (N+2)*ny)

cu.csSolve!(wk, m, opt, traj0, params0, cu.WRT_TRAJ)

# # setup opt ---

# cu.Jobj(m, traj0, params0)
# # DJ = similar(traj0)
# # cu.∇Jobj!(DJ, m, traj0, params0)

# # gL, gU = cu.gbounds(m, traj0)
# # g0 = similar(gL)
# # cu.gvalues!(g0, m, traj0, params0)
# # nnz = cu.Dgnnz(m, N) # 532
# # row = zeros(Int32, nnz)
# # col = similar(row)
# # val = zeros(nnz)
# # # cu.Dgsparse!(row, col, val, m, traj0, params0, true)
# # cu.Dgsparse!(row, col, val, m, traj0, params0, :Structure)

# # IPOPT
# # prob = cu.nloptsetup(m, traj0, params0; fixedδt=0.3)
# # status = cu.nloptsolve(prob)
# # Ipopt.ApplicationReturnStatus[status]

# trajs, params = cu.csAlternateSolve(m, opt, traj0, params0, 2; μst=[1e-2,1e-2], Ninnert=2, μsp=[1e-2,1e-2], Ninnerp=2)

# pl1 = plotTrajs(m, opt, trajt, params0, (trajs[:,i] for i = 1:size(trajs,2))...)
# pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# display(params)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))
# gui()
