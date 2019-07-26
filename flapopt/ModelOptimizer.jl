
#=========================================================================
Dynamics constraint
=========================================================================#

"""
Dynamics constraint at state y, input u
	
Note that the actual constraints are:
	g[1:ny] = y1 - y1(0) initial condition
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
		gout[vy2] = traj[vy2] - (traj[vy] + δt * dydt(m, traj[vy], traj[vu], params))
	end

	# Initial condition
	gout[liy[:,1]] = -traj[@view liy[:,1]]

	return
end

# "Bounds corresponding to the constraints above"
function gbounds(m::Model, N::Int)::Tuple{Vector, Vector}
	ny, nu = dims(m)
    g_L = zeros((N+1)*ny)
    g_U = similar(g_L)
    # first N*ny = 0 (dynamics)
    return g_L, g_U
end

function Dgnnz(m::Model, N::Int; vart::Bool=true)::Int
	ny, nu = dims(m)
	# Assuming the Jacobians are dense. The terms below correspond to the "-I"s, the "I + δt A"'s, the "δt B"'s
	return ny*(N+1) + ny^2*N + ny*nu*N
end

function Dgsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, traj::Vector, params::Vector, setVals::Bool; vart::Bool=true, fixedδt::Float64=1e-3, order::Int=1)
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	liy, liu = linind(m, N)
	δt = vart ? traj[end] : fixedδt
	# Preallocate outputs
	df_dy = zeros(ny, ny)
	df_du = zeros(ny, nu)

	# Fill in -I's
	for ii = 1:ny*(N+1)
		if setVals
			value[ii] = -1
		else
			row[ii] = col[ii] = ii;
		end
	end
	
	# Offsets into the row[], col[], val[] arrays. These will be incremented and keep track of the index in the loops below.
	offsA = ny*(N+1)
	offsB = ny*(N+1) + ny^2*N

	# Fill in Jacobians
	for k = 1:N
		# Get the jacobians at this y, u
		Df!(df_dy, df_du, m, traj[@view liy[:,k]], traj[@view liu[:,k]], params)

		# Insert A NOTE j outer loop for Julia's col-major storage and better loop unrolling
		for j = 1:ny
			for i = 1:ny
				if setVals
					value[offsA] = δt * df_dy[i,j] + (i == j ? 1 : 0)
				else
					row[offsA] = k*ny + i
					col[offsA] = (k-1)*ny + j
				end
				offsA += 1
			end
		end

		# Insert B
		for j = 1:nu
			for i = 1:ny
				if setVals
					value[offsB] = δt * df_du[i,j]
				else
					row[offsA] = k*ny + i
					col[offsA] = (N+k)*ny + j
				end
				offsB += 1
			end
		end
	end
	return
end

#=========================================================================
Objective
=========================================================================#

function ∇Jobj!(∇Jout, m::Model, traj::Vector, params::Vector; vart::Bool=true)
	Jtraj(tt::Vector) = Jobj(m, tt, params; vart=vart)
	ForwardDiff.gradient!(∇Jout, Jtraj, traj)
	return
end

#=========================================================================
Solver interface
=========================================================================#

function nloptsetup(m::Model, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3)
	ny, nu = dims(m)
	# Construct constraints
	N = Nknot(m, traj; vart=vart)

	# Define the things needed for IPOPT
	x_L, x_U = xbounds(m, N; vart=vart)
	g_L, g_U = gbounds(m, N)
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, x, params; vart=vart, fixedδt=fixedδt)
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = Dgsparse!(rows, cols, values, m, x, params, mode == :Values; vart=vart, fixedδt=fixedδt)
	eval_f(x::Vector{Float64}) = Jobj(m, x, params; vart=vart, fixedδt=fixedδt)
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = ∇Jobj!(grad_f, x)

	# Create IPOPT problem
	prob = Ipopt.createProblem(
		length(traj), # Number of variables
		x_L, # Variable lower bounds
		x_U, # Variable upper bounds
		length(g_L), # Number of constraints
		g_L,       # Constraint lower bounds
		g_U,       # Constraint upper bounds
		Dgnnz(m, N; vart=vart),  # Number of non-zeros in Jacobian
		0,             # Number of non-zeros in Hessian
		eval_f,                     # Callback: objective function
		eval_g,                     # Callback: constraint evaluation
		eval_grad_f,                # Callback: objective function gradient
		eval_jac_g,                 # Callback: Jacobian evaluation
		nothing           # Callback: Hessian evaluation
	)
	Ipopt.addOption(prob, "hessian_approximation", "limited-memory")

	return prob
end

