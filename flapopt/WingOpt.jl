
using ForwardDiff, DifferentialEquations, Plots
include("Wing2DOF.jl")

y0 = [1e-2, 0, 0, 0]
u0 = [0.]
params0 = [0.005,1.5]
paeroFun(q::Vector) = Wing2DOF.aero([q;[0,0]], u0, params0)[1]
println("paero ", paeroFun(y0[1:2]))

# Try to use ForwardDiff
JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# println(JaeroFun(y0[1:2]))
println(Wing2DOF.aero(y0, u0, params0)[2])

# Create a traj
function strokePosControlVF(y, p, t)
	σdes = 0.75 * 0.6 * sin(150 * 2 * π * t)
	τ = 1e2 * (σdes - y[1]) - 1e-2 * y[3]
	return Wing2DOF.dydt(y, [τ], params0)
end

Y = [0.963506  0.00615112  0.00404816  0.309899  0.466694;
 0.417588  0.452104    0.0708825   0.269745  0.796851;
 0.808803  0.0559587   0.0637219   0.174673  0.607593;
 0.548119  0.397353    0.61362     0.538941  0.166001]
for i=1:5
	println(strokePosControlVF(Y[:,i], [], 0))
end

# println(strokePosControlVF(y0[3:4], y0[1:2], [], 0))


# dy = zeros(4)
# println(dy)
# wing2dof!(dy, y0, [], 0)
# println(dy)

# tspan = (0.0, 0.01)
# prob = ODEProblem{false}(strokePosControlVF, y0, tspan)
# sol = solve(prob)

# println("SOLVED!")

# plotly()
# # # println(sol)
# plot(sol[1:2,:])

# gui()
