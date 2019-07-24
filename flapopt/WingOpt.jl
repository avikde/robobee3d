
using DifferentialEquations, Plots
gr() # backend
include("Wing2DOF.jl")
include("WingOptimizer.jl")

y0 = [1e-2, 0, 0, 0]
u0 = [0.]
params0 = [0.005,1.5]

# paeroFun(q::Vector) = Wing2DOF.aero([q;[0,0]], u0, params0)[1]
# println("paero ", paeroFun(y0[1:2]))

# # Try to use ForwardDiff
# JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# # println(JaeroFun(y0[1:2]))
# println(Wing2DOF.aero(y0, u0, params0)[2])

# Create a traj
function strokePosController(y, t)
	σdes = 0.75 * 0.6 * sin(150 * 2 * π * t)
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
olTrajt = sol.t[olRange]
N = length(olTrajt)
olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
olTraju = [strokePosController(olTrajaa[i], olTrajt[i]) for i in 1:N-1] # get u1,...,u(N-1)
traj0 = [vcat(olTrajaa...); olTraju] # dirtran form {x1,..,xN,u1,...,u(N-1)}

