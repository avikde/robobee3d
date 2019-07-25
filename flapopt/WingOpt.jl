
using Plots, BenchmarkTools, StaticArrays
gr() # backend
include("Wing2DOF.jl")
include("WingOptimizer.jl")

y0 = [1e-2, 0, 0, 0]
u0 = [0.]
params0 = [0.005, 1.5]

# paeroFun(q::Vector) = Wing2DOF.aero([q;[0,0]], u0, params0)[1]
# println("paero ", paeroFun(y0[1:2]))

# # Try to use ForwardDiff
# JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# # println(JaeroFun(y0[1:2]))
# println(Wing2DOF.aero(y0, u0, params0)[2])

trajt, traj0 = Wing2DOF.createInitialTraj(150, [1e2, 1e-2], params0)

WingOptimizer.plotTrajs(trajt, traj0)

ly, lu = WingOptimizer.linind(traj0)
# println(traj0[ly[:,2]], traj0[lu[:,2]])
paero = @SVector zeros(2)
Faero = @SVector zeros(2)
Jaero = @SMatrix zeros(2,2)
# @btime Wing2DOF.aero!(paero, Jaero, Faero, y0, u0, params0)
# @btime Wing2DOF.eval_f(traj0, params0)
# println()

g0 = zeros(ly[WingOptimizer.ny,N])
@btime Wing2DOF.eval_g!(traj0, g0, params0)
# println(g0)
