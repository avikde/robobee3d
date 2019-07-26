
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using Plots, BenchmarkTools
gr() # backend
using Revise # while developing
import controlutils
cu = controlutils
include("Wing2DOF.jl")

# create an instance
m = Wing2DOFModel()
ny, nu = cu.dims(m)
N = 23
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

# cu.Jobj(m, traj0, params0)
DJ = similar(traj0)
@btime cu.∇Jobj!(DJ, m, traj0, params0)

# # @btime cu.Df(m, [0.1,0,0,0], [10.], params0)
# g0 = zeros(N*ny)
# # cu.g_LU(m, N)
# @btime cu.eval_g!(g0, m, traj0, params0)
# row = Int32[]
# col = Int32[]
# vals = Float64[]
# @btime cu.Dgsparse!(row, col, vals, m, traj0, params0, true)
# cu.nloptsetup(m, traj0, params0)
# # @btime Wing2DOF.eval_g!(traj0, params0, N, ly, lu, gout)
# # println(g0)

