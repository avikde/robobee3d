
include("Model.jl") #< including this helps vscode reference the functions in there

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
function gbounds(m::Model, traj::Vector; vart::Bool=true)::Tuple{Vector, Vector}
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	# println("CALLED gbounds with $(ny) $(nu) $(N) $(pointer_from_objref(traj))")
	g_L = [-traj[1:ny]; zeros(N*ny)]
    g_U = [-traj[1:ny]; zeros(N*ny)]
    # first N*ny = 0 (dynamics)
    return g_L, g_U
end

function Dgnnz(m::Model, N::Int; vart::Bool=true)::Int
	ny, nu = dims(m)
	# Assuming the Jacobians are dense. The terms below correspond to the "-I"s, the "I + δt A"'s, the "δt B"'s
	return ny*(N+1) + ny^2*N + ny*nu*N
end

function Dgsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, traj::Vector, params::Vector, mode; vart::Bool=true, fixedδt::Float64=1e-3, order::Int=1)
	ny, nu = dims(m)
	N = Nknot(m, traj; vart=vart)
	liy, liu = linind(m, N)

	# NOTE: traj is NULL when in :Structure mode
	if mode != :Structure
		δt = vart ? traj[end] : fixedδt
		# Preallocate outputs
		df_dy = zeros(ny, ny)
		df_du = zeros(ny, nu)
	end

	# Fill in -I's
	for ii = 1:ny*(N+1)
		if mode == :Structure
			row[ii] = ii
			col[ii] = ii
		else
			value[ii] = -1
		end
	end
	
	# Offsets into the row[], col[], val[] arrays. These will be incremented and keep track of the index in the loops below.
	offsA = ny*(N+1)
	offsB = ny*(N+1) + ny^2*N

	# Fill in Jacobians
	for k = 1:N
		if mode != :Structure
			# Get the jacobians at this y, u
			Df!(df_dy, df_du, m, traj[@view liy[:,k]], traj[@view liu[:,k]], params)
		end

		# Insert A
		# NOTE: j outer loop for Julia's col-major storage and better loop unrolling
		for j = 1:ny
			for i = 1:ny
				offsA += 1 # needs to be up here due to 1-indexing!
				if mode == :Structure
					row[offsA] = k*ny + i
					col[offsA] = (k-1)*ny + j
				else
					value[offsA] = δt * df_dy[i,j] + (i == j ? 1 : 0)
				end
			end
		end

		# Insert B
		for j = 1:nu
			for i = 1:ny
				offsB += 1
				if mode == :Structure
					row[offsB] = k*ny + i
					col[offsB] = (N+1)*ny + (k-1)*nu + j
				else
					value[offsB] = δt * df_du[i,j]
				end
			end
		end
	end
	return
end

#=========================================================================
Objective
=========================================================================#

function ∇Jobj!(∇Jout::Vector, m::Model, traj::Vector, params::Vector; vart::Bool=true)
	# println("CALLED DJobj with $(pointer_from_objref(traj))")
	Jtraj(tt::Vector) = Jobj(m, tt, params; vart=vart)
	ForwardDiff.gradient!(∇Jout, Jtraj, traj)
	return
end

#=========================================================================
Custom solver
=========================================================================#

function mysol(m::Model, traj0::Vector, params::Vector; μs::Array{Float64}=[1e-1], Ninner::Int=1, vart::Bool=true, fixedδt::Float64=1e-3)
	ny, nu = dims(m)
	N = Nknot(m, traj0; vart=vart)
	liy, liu = linind(m, N)
	δt = vart ? traj0[end] : fixedδt

	# Constraint bounds
	x_L, x_U = xbounds(m, N; vart=vart)
	g_L, g_U = gbounds(m, traj0; vart=vart)

	# Preallocate outputs
	g = similar(g_L)
	Nx, Ng = length(traj0), length(g)
	df_dy = zeros(ny, ny)
	df_du = zeros(ny, nu)
	Ak = similar(df_dy)
	Bk = similar(df_du)
	∇J = zeros(Nx)

	# TODO: sparse matrices for these spzeros
	HJ = zeros(Nx, Nx)
	DgTg = zeros(Nx) # Dg' * g
	# TODO: better Dg' Dg computation that doesn't compute Dg
	Dg = zeros(Ng, Nx)
	for ii = 1:(N+1)*ny
		Dg[ii,ii] = -1
	end
	
	# Take some number of steps
	traj = copy(traj0)
	traj1 = similar(traj)
	for μ in μs
		for stepi = 1:Ninner
			println("μ=$(μ), step=$(stepi)")
			# One step
			gvalues!(g, m, traj, params; vart=vart, fixedδt=fixedδt)
			for k = 1:N
				Df!(df_dy, df_du, m, traj[@view liy[:,k]], traj[@view liu[:,k]], params)
				# Dg_y' * g = [-g0 + A1^T g1, ..., -g(N-1) + AN^T gN, -gN]
				Ak .= I + δt * df_dy
				DgTg[liy[:,k]] .= -g[@view liy[:,k]] + Ak' * g[@view liy[:,k+1]]
				# Dg_u' * g = [B1^T g1, ..., BN^T gN]
				Bk .= δt * df_du
				DgTg[liu[:,k]] .= Bk' * g[@view liy[:,k+1]]
				# Dg_δt' * g = 0

				# TODO: better Dg' Dg computation that doesn't compute Dg
				Dg[liy[:,k+1], liy[:,k]] .= Ak
				Dg[liy[:,k+1], liu[:,k]] .= Bk
			end
			DgTg[liy[:,N+1]] .= -g[@view liy[:,N+1]]
			# Calculate cost gradient from objective and an added penalty term
			∇Jobj!(∇J, m, traj, params; vart=vart)
			∇J .= ∇J + μ * DgTg

			# Hessian of objective and add penalty term
			ForwardDiff.hessian!(HJ, tt::Vector -> Jobj(m, tt, params; vart=vart), traj)
			HJ .= 1/2 * (HJ' + HJ) # take the symmetric part

			# TODO: better Dg' Dg computation that doesn't compute Dg
			DgTDg = Dg' * Dg
			HJ .= HJ + μ * DgTDg

			# # Gradient descent
			# v = -∇J
			# Newton or Gauss-Newton. Use PositiveFactorizations.jl to ensure psd Hessian
			v = -(cholesky(Positive, HJ) \ ∇J)

			function Jcallable(x::Vector)::Float64
				gvalues!(g, m, x, params; vart=vart, fixedδt=fixedδt)
				return Jobj(m, x, params; vart=vart, fixedδt=fixedδt) + μ/2 * g' * g
			end

			_backtrackingLineSearch!(traj1, traj, ∇J, v, Jcallable; α=0.2, β=0.7)
			traj .= traj1
		end
	end
	return traj
end

function _backtrackingLineSearch!(x1::Vector, x0::Vector, ∇J0::Vector, v::Vector, Jcallable; α::Float64=0.45, β::Float64=0.9)
	σ = 1
	J0 = Jcallable(x0)
	# search for step size
	while true
		σ = β * σ
		x1 .= x0 + σ * v
		J1 = Jcallable(x1)
		# debug line search
		println("J0=$(round(J0; sigdigits=4)), J1=$(round(J1; sigdigits=4)), σ=$(round(σ; sigdigits=6))")
		if J1 < J0 + α * σ * ∇J0' * v || σ < 1e-6
			return J1
		end
	end
	return J0
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
	g_L, g_U = gbounds(m, traj; vart=vart)
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, x, params; vart=vart, fixedδt=fixedδt)
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = Dgsparse!(rows, cols, values, m, x, params, mode; vart=vart, fixedδt=fixedδt)
	eval_f(x::Vector{Float64}) = Jobj(m, x, params; vart=vart, fixedδt=fixedδt)
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = ∇Jobj!(grad_f, m, x, params)

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

	# TODO: this should be an update only without need to setup
	prob.x = traj

	return prob
end

nloptsolve(prob) = Ipopt.solveProblem(prob)
