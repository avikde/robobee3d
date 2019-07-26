
using Ipopt

include("Model.jl")

#===========================================================================
Dynamics constraint
===========================================================================#

"""
Dynamics constraint at state y, input u
	
Note that the actual constraints are:
	g = -ynext + (y + δt * fy)
	dg_dynext = -I
	dg_dy = δt * df_dy + I
	dg_du = δt * df_du
	dg_dδt = fy
"""
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

function DgsparseNNZ(m::Model, N::Int)::Int
	# TODO:
	return 1
end

function Dgsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, traj::Vector, params::Vector, setVals::Bool; vart::Bool=true, order::Int=1)
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	liy, liu = linind(m, N)
	# Preallocate outputs
	df_dy = zeros(ny, ny)
	df_du = zeros(ny, nu)

	for k = 1:N
        vy = @view liy[:,k]
		vu = @view liu[:,k]
		df_dy[:], df_du[:] = Df(m, traj[vy], traj[vu], params)
		# TODO:
		if setVals

		else

		end
	end
end

#===========================================================================
Objective
===========================================================================#

function DJobj!(DJ::Vector, traj::Vector)
	# FIXME:
	DJ = similar(traj)
end

#===========================================================================
Solver interface
===========================================================================#

function nloptsetup(m::Model, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3)
	ny, nu = dims(m)
	# Construct constraints
	N = Nknot(m, traj; vart=vart)

	# Define the things needed for IPOPT
	x_L, x_U = xbounds(m, N; vart=vart)
	g_L, g_U = gbounds(m, N)
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, x, params; vart=vart, fixedδt=fixedδt)
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = Dgsparse!(rows, cols, values, m, x, params, mode == :Values; vart=vart)
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = DJobj!(grad_f, x)

	# Create IPOPT problem
	prob = createProblem(
		length(traj), # Number of variables
		x_L, # Variable lower bounds
		x_U, # Variable upper bounds
		length(g_L), # Number of constraints
		g_L,       # Constraint lower bounds
		g_U,       # Constraint upper bounds
		DgsparseNNZ(m, N),  # Number of non-zeros in Jacobian
		0,             # Number of non-zeros in Hessian
		eval_f,                     # Callback: objective function
		eval_g,                     # Callback: constraint evaluation
		eval_grad_f,                # Callback: objective function gradient
		eval_jac_g,                 # Callback: Jacobian evaluation
		nothing           # Callback: Hessian evaluation
	)
	addOption(prob, "hessian_approximation", "limited-memory")

	return prob
end

