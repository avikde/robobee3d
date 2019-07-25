
using DifferentialEquations, Plots, BenchmarkTools, StaticArrays
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

# Create a traj
σmax = Wing2DOF.limits()[end][1]
function strokePosController(y, t)
	σdes = 0.75 * σmax * sin(150 * 2 * π * t)
	return 1e2 * (σdes - y[1]) - 1e-2 * y[3]
end
strokePosControlVF(y, p, t) = Wing2DOF.dydt(y, [strokePosController(y, t)], params0)
# OL traj1
teval = collect(0:1e-4:0.1)
prob = ODEProblem(strokePosControlVF, y0, (teval[1], teval[end]))
sol = solve(prob, saveat=teval)

# σt = plot(sol, vars=1, linewidth=1, ylabel="act disp [m]")
# Ψt = plot(sol, vars=2, linewidth=1, ylabel="hinge ang [r]")
# plot(σt, Ψt, layout=(2,1))
# gui()

olRange = 171:3:238
trajt = sol.t[olRange]
δt = trajt[2] - trajt[1]
N = length(trajt) - 1
olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
olTraju = [strokePosController(olTrajaa[i], trajt[i]) for i in 1:N] # get u1,...,uN
traj0 = [vcat(olTrajaa...); olTraju; δt] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

WingOptimizer.plotTrajs(trajt, traj0)

# ly, lu = WingOptimizer.linind(traj0)
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
