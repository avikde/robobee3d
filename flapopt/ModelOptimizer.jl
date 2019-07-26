
using Ipopt

include("Model.jl")

function gvalues!(gout::Vector, m::Model, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3, order::Int=1)
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	liy, liu = linind(m, N)
	δt = vart ? traj[end] : fixedδt

	for k = 1:N
        vy2 = @view liy[:,k+1]
        vy = @view liy[:,k]
		vu = @view liu[:,k]
		gout[vy] = traj[vy2] - (traj[vy] + δt * dydt(m, traj[vy], traj[vu], params))
	end
end

# "Bounds corresponding to the constraints above"
function gbounds(m::Model, N::Int)::Tuple{Vector, Vector}
	ny, nu = dims(m)
    g_L = zeros(N*ny)
    g_U = similar(g_L)
    # first N*ny = 0 (dynamics)
    return g_L, g_U
end


function eval_jac_g!(m::Model, traj::Vector, params::Vector; vart::Bool=true, order::Int=1)::Tuple
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	liy, liu = linind(m, N)
	# Preallocate outputs
	fy = zeros(ny)
	# df_dy = zeros(ny, ny)
	# df_du = zeros(ny, nu)

	for k = 1:N
        vy2 = @view liy[:,k+1]
        vy = @view liy[:,k]
		vu = @view liu[:,k]
		g[vy] = traj[vy2] - (traj[vy] + δt * dydt(traj[vy], traj[vu], params))
		# Get the matrices
		cu.gdyn!(fy, df_dy, df_du, m, traj[vy], traj[vu], params0)
	end
end


function nloptsetup(m::Model, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3)
	ny, nu = dims(m)
	# Construct constraints
	N = Nknot(m, traj; vart=vart)

	# Define the things needed for IPOPT
	x_L, x_U = xbounds(m, N; vart=vart)
	g_L, g_U = gbounds(m, N)
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, x, params; vart=vart, fixedδt=fixedδt)

	# prob = createProblem(
	# 	length(traj), # Number of variables
	# 	x_L, # Variable lower bounds
	# 	x_U, # Variable upper bounds
	# 	length(g_L), # Number of constraints
	# 	g_L,       # Constraint lower bounds
	# 	g_U,       # Constraint upper bounds
	# 	nele_jac::Int,              # Number of non-zeros in Jacobian
	# 	0,             # Number of non-zeros in Hessian
	# 	eval_f,                     # Callback: objective function
	# 	eval_g!,                     # Callback: constraint evaluation
	# 	eval_grad_f,                # Callback: objective function gradient
	# 	eval_jac_g,                 # Callback: Jacobian evaluation
	# 	eval_h = nothing           # Callback: Hessian evaluation
	# )
	# addOption(prob, "hessian_approximation", "limited-memory")

	# return prob
end

