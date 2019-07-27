
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

gL, gU = cu.gbounds(m, traj0)
g0 = similar(gL)
cu.gvalues!(g0, m, traj0, params0)
nnz = cu.Dgnnz(m, N) # 532
row = zeros(Int32, nnz)
col = similar(row)
val = zeros(nnz)
# cu.Dgsparse!(row, col, val, m, traj0, params0, true)
cu.Dgsparse!(row, col, val, m, traj0, params0, :Structure)

println(pointer_from_objref(traj0))

# println(pointer_from_objref(row), pointer_from_objref(col), pointer_from_objref(val))
# eval_jac_g(traj0, :Structure, row, col, val)
# println(pointer_from_objref(row), pointer_from_objref(col), pointer_from_objref(val))



# Define the things needed for IPOPT
x_L, x_U = cu.xbounds(m, N)
g_L, g_U = cu.gbounds(m, traj0)
eval_g(x::Vector, g::Vector) = cu.gvalues!(g, m, x, params0)
eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = cu.Dgsparse!(rows, cols, values, m, x, params0, mode)
eval_f(x::Vector{Float64}) = cu.Jobj(m, x, params0)
eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = cu.∇Jobj!(grad_f, m, x, params0)

eval_jac_g(traj0, :Structure, row, col, val)

# # Create IPOPT problem
# prob = Ipopt.createProblem(
# 	length(traj0), # Number of variables
# 	x_L, # Variable lower bounds
# 	x_U, # Variable upper bounds
# 	length(g_L), # Number of constraints
# 	g_L,       # Constraint lower bounds
# 	g_U,       # Constraint upper bounds
# 	cu.Dgnnz(m, N),  # Number of non-zeros in Jacobian
# 	0,             # Number of non-zeros in Hessian
# 	eval_f,                     # Callback: objective function
# 	eval_g,                     # Callback: constraint evaluation
# 	eval_grad_f,                # Callback: objective function gradient
# 	eval_jac_g,                 # Callback: Jacobian evaluation
# 	nothing           # Callback: Hessian evaluation
# )
# prob.x = traj0

# Ipopt.addOption(prob, "hessian_approximation", "limited-memory")

# println(pointer_from_objref(prob.x))

# # Ipopt.solveProblem(prob)


# prob = cu.nloptsetup(m, traj0, params0; fixedδt=0.3)
# status = cu.nloptsolve(prob)
# Ipopt.ApplicationReturnStatus[status]
