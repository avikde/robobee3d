
# This is a "script" to test/run the functions from

using Plots, BenchmarkTools, StaticArrays
gr() # backend
include("Wing2DOF.jl")
include("ModelOptimizer.jl")

# create an instance
m = Wing2DOFModel()
ny, nu = dims(m)
N = 23
liy, liu = linind(m, N)

params0 = [2, 20] # cbar, T

# println("paero ", paeroFun(y0[1:2]))

# # Try to use ForwardDiff
# using ForwardDiff
# paeroFun(q::Vector) = w2daero([q;[0,0]], [10], params0)[1]
# @btime JaeroFun = q -> ForwardDiff.jacobian(paeroFun, q)
# @btime J1 = JaeroFun([0.1,1])
# # println()
# println(Wing2DOF.aero(y0, u0, params0)[2])

trajt, traj0 = createInitialTraj(m, N, 0.15, [1e3, 1e2], params0)
plotTrajs(m, trajt, params0, traj0)

# println(traj0[ly[:,2]], traj0[lu[:,2]])
paero = @SVector zeros(2)
Faero = @SVector zeros(2)
Jaero = @SMatrix zeros(2,2)
# @btime Wing2DOF.aero!(paero, Jaero, Faero, y0, u0, params0)
# @btime Wing2DOF.eval_f(traj0, params0)
# println()

# setup opt ---

optsetup(m, traj0)
# @btime Wing2DOF.eval_g!(traj0, params0, N, ly, lu, gout)
# println(g0)
