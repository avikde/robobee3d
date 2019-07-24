
using ForwardDiff, DifferentialEquations, Plots
include("Wing2DOF.jl")

y0 = [1e-2, 0, 0, 0]
u0 = [0.]
params0 = [0.05,1.5]
paeroFun(q::Vector) = Wing2DOF.aero([q;[0,0]], u0, params0)[1]
println("paero ", paeroFun(y0[1:2]))

# Try to use ForwardDiff
JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# println(JaeroFun(y0[1:2]))
println(Wing2DOF.aero(y0, u0, params0)[2])

# Create a traj
function strokePosControlVF(y, p, t)
	σdes = 0.75 * 0.15 * sin(150 * 2 * π * t)
	τ = 1e2 * (σdes - y[1]) - 1e-2 * y[3]
	return Wing2DOF.dydt(y, [τ], params0)
end

# println(strokePosControlVF(y0[3:4], y0[1:2], [], 0))


# dy = zeros(4)
# println(dy)
# wing2dof!(dy, y0, [], 0)
# println(dy)

tspan = (0.0, 0.01)
prob = ODEProblem{false}(strokePosControlVF, y0, tspan)
sol = solve(prob)

println("SOLVED!")

# # println(sol)
plot(sol[1:2,:])

gui()
