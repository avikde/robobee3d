
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using Plots, BenchmarkTools
import Ipopt
gr() # backend
using Revise # while developing
import controlutils
cu = controlutils
include("Wing2DOF.jl")

# create an instance
m = Wing2DOFModel()
ny, nu = cu.dims(m)
N = 22
liy, liu = cu.linind(m, N)

params0 = [2.0, 20.0] # cbar, T

# println("paero ", paeroFun(y0[1:2]))

# # Try to use ForwardDiff
# using ForwardDiff
# paeroFun(q::Vector) = w2daero([q;[0,0]], [10], params0)[1]
# @btime JaeroFun = q -> ForwardDiff.jacobian(paeroFun, q)
# @btime J1 = JaeroFun([0.1,1])
# # println()
# println(Wing2DOF.aero(y0, u0, params0)[2])

trajt, traj0 = createInitialTraj(m, N, 0.15, [1e3, 1e2], params0)
# cu.plotTrajs(m, trajt, params0, traj0)

# setup opt ---

cu.Jobj(m, traj0, params0)
# DJ = similar(traj0)
# cu.∇Jobj!(DJ, m, traj0, params0)

# gL, gU = cu.gbounds(m, traj0)
# g0 = similar(gL)
# cu.gvalues!(g0, m, traj0, params0)
# nnz = cu.Dgnnz(m, N) # 532
# row = zeros(Int32, nnz)
# col = similar(row)
# val = zeros(nnz)
# # cu.Dgsparse!(row, col, val, m, traj0, params0, true)
# cu.Dgsparse!(row, col, val, m, traj0, params0, :Structure)

prob = cu.nloptsetup(m, traj0, params0; fixedδt=0.3)
status = cu.nloptsolve(prob)
Ipopt.ApplicationReturnStatus[status]
