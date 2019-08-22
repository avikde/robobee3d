
include("OptBase.jl") #< including this helps vscode reference the functions in there

#=========================================================================
IPOPT Solver interface
=========================================================================#

# "Bounds corresponding to the constraints above"
function gbounds(m::Model, opt::OptOptions, traj::Vector, εic::Float64=0., εdyn::Float64=0., εsymm::Float64=0.)::Tuple{Vector, Vector}
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)

	# println("CALLED gbounds with $(ny) $(nu) $(N) $(pointer_from_objref(traj))")
	gU = [fill(εic, ny); fill(εdyn, N*ny)]
	if opt.boundaryConstraint == :symmetric
		gU = [gU; fill(εsymm, ny)]
	end
    return -gU, gU
end

function Dgnnz(m::Model, opt::OptOptions, traj::Vector)::Int
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# Assuming the Jacobians are dense. The terms below correspond to initial cond., dynamics
	nnz = ny + (ny + ny^2)*N + ny*nu*N
	if opt.boundaryConstraint == :symmetric
		nnz += 2*ny # the two I's for symmetry
	end
	return nnz
end

function Dgsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, opt::OptOptions, traj::Vector, params::Vector, mode, ny, nu, N, δt)
	#=
	NOTE ON USAGE

	If the iRow and jCol arguments are not NULL, then IPOPT wants you to fill in the sparsity structure of the Jacobian (the row and column indices only). At this time, the x argument and the values argument will be NULL.
	If the x argument and the values argument are not NULL, then IPOPT wants you to fill in the values of the Jacobian as calculated from the array x (using the same order as you used when specifying the sparsity structure). At this time, the iRow and jCol arguments will be NULL;
	=#
	# ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	if mode != :Structure
		ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
		yk = k -> @view traj[liy[:,k]]
		uk = k -> @view traj[liu[:,k]]
		# Preallocate outputs
		Ak = zeros(ny, ny)
		Bk = zeros(ny, nu)
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
			dlin!(Ak, Bk, m, yk(k), uk(k), params, δt)
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
					value[offsA] = Ak[i,j]
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
					value[offsB] = Bk[i,j]
				end
			end
		end
	end

	# the new symmetry constraint
	if opt.boundaryConstraint == :symmetric
		for j = 1:ny
			# Two -I's
			offsB += 1
			if mode == :Structure
				row[offsB] = (N+1) * ny + j
				col[offsB] = j
			else
				value[offsB] = -1
			end
			offsB += 1
			if mode == :Structure
				row[offsB] = (N+1) * ny + j
				col[offsB] = (N) * ny + j
			else
				value[offsB] = -1
			end
		end
	end
	return
end

#=========================================================================
Param opt Dg
=========================================================================#

function Dgpsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, opt::OptOptions, traj::Vector, params::Vector, mode, ny, nu, N, δt)
	np = pdims(m)
	if mode != :Structure
		ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
		yk = k -> @view traj[liy[:,k]]
		uk = k -> @view traj[liu[:,k]]
		# Preallocate outputs
		Pk = zeros(ny, np)
	end
	offs = 0
	
	# Fill in Jacobians
	for k = 1:N
		if mode != :Structure
			# Get the jacobians at this y, u
			dlinp!(Pk, m, yk(k), uk(k), params, δt)
		end

		# Insert A
		# NOTE: j outer loop for Julia's col-major storage and better loop unrolling
		for j = 1:np
			for i = 1:ny
				offs += 1 # needs to be up here due to 1-indexing!
				if mode == :Structure
					row[offs] = k*ny + i
					col[offs] = j
				else
					value[offs] = Pk[i,j]
				end
			end
		end
	end
	return
end

#=========================================================================
Direct Collocation
=========================================================================#

function DgnnzDirCol(m::Model, opt::OptOptions, traj::Vector)::Int
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# Assuming the Jacobians are dense.
	nnz = ny + ny^2*2*N + ny*nu*2*N
	if opt.boundaryConstraint == :symmetric
		nnz += 2*ny # the two I's for symmetry
	end
	return nnz
end

function DgsparseDirCol!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, opt::OptOptions, traj::Vector, params::Vector, mode, ny, nu, N, δt)
	#=
	NOTE ON USAGE

	If the iRow and jCol arguments are not NULL, then IPOPT wants you to fill in the sparsity structure of the Jacobian (the row and column indices only). At this time, the x argument and the values argument will be NULL.
	If the x argument and the values argument are not NULL, then IPOPT wants you to fill in the values of the Jacobian as calculated from the array x (using the same order as you used when specifying the sparsity structure). At this time, the iRow and jCol arguments will be NULL;
	=#
	# ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	if mode != :Structure
		ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
		yk = k -> @view traj[liy[:,k]]
		uk = k -> @view traj[liu[:,k]]
		# Preallocate outputs
		dg_dyk = zeros(ny, ny)
		dg_dykp1 = zeros(ny, ny)
		dg_duk = zeros(ny, nu)
		dg_dukp1 = zeros(ny, nu)
	end

	# For IC
	for ii = 1:ny
		if mode == :Structure
			row[ii] = ii
			col[ii] = ii
		else
			value[ii] = -1
		end
	end
	
	# Offsets into the row[], col[], val[] arrays. These will be incremented and keep track of the index in the loops below.
	offs = ny

	# Fill in Jacobians
	for k = 1:N
		if mode != :Structure
			if k < N
				ukp1 = uk(k+1)
			elseif k == N && opt.boundaryConstraint == :symmetric
				ukp1 = -uk(1)
			else
				ukp1 = uk(k) # don't have any more u's in the current parameterization
			end
			# Get the jacobians at this y, u
			ForwardDiff.jacobian!(dg_dyk, yy -> gkDirCol(m, params, yy, yk(k+1), uk(k), ukp1, δt)[1], yk(k))
			ForwardDiff.jacobian!(dg_dykp1, yy -> gkDirCol(m, params, yk(k), yy, uk(k), ukp1, δt)[1], yk(k+1))
			ForwardDiff.jacobian!(dg_duk, uu -> gkDirCol(m, params, yk(k), yk(k+1), uu, ukp1, δt)[1], uk(k))
			ForwardDiff.jacobian!(dg_dukp1, uu -> gkDirCol(m, params, yk(k), yk(k+1), uk(k), uu, δt)[1], ukp1)
		end

		# Insert A
		# NOTE: j outer loop for Julia's col-major storage and better loop unrolling
		for j = 1:ny
			for i = 1:ny
				# dgk/dyk
				offs += 1 # needs to be up here due to 1-indexing!
				if mode == :Structure
					row[offs] = k*ny + i
					col[offs] = (k-1)*ny + j
				else
					value[offs] = dg_dyk[i,j]
				end
				# dgk/dykp1
				offs += 1 # needs to be up here due to 1-indexing!
				if mode == :Structure
					row[offs] = k*ny + i
					col[offs] = k*ny + j
				else
					value[offs] = dg_dykp1[i,j]
				end
			end
		end

		# Insert B
		for j = 1:nu
			for i = 1:ny
				# dgk/duk
				offs += 1
				if mode == :Structure
					row[offs] = k*ny + i
					col[offs] = (N+1)*ny + (k-1)*nu + j
				else
					value[offs] = dg_duk[i,j]
				end
				# dgk/dukp1
				offs += 1
				if mode == :Structure
					row[offs] = k*ny + i
					# FIXME: will go out of range
					col[offs] = (N+1)*ny + k*nu + j
				else
					value[offs] = dg_dukp1[i,j]
				end
			end
		end
	end

	# the new symmetry constraint
	if opt.boundaryConstraint == :symmetric
		for j = 1:ny
			# Two -I's
			offs += 1
			if mode == :Structure
				row[offs] = (N+1) * ny + j
				col[offs] = j
			else
				value[offs] = -1
			end
			offs += 1
			if mode == :Structure
				row[offs] = (N+1) * ny + j
				col[offs] = (N) * ny + j
			else
				value[offs] = -1
			end
		end
	end
	return
end


"""
IPOPT problem struct (output) fields:
Solution of the constraint multipliers, lambda = prob.mult_g
Solution of the bound multipliers, z_L and z_U = prob.mult_x_L, prob.mult_x_U
Objective value f(x*) = prob.obj_val
Final constraint values = prob.g
"""
function ipoptsolve(m::Model, opt::OptOptions, traj::Vector, params::Vector, εs, optWrt::Symbol; kwargs...)
	optWrt in OptVar || throw(ArgumentError("invalid optWrt: $optWrt"))

	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	np = pdims(m)
	δt = opt.vart ? traj[end] : opt.fixedδt

	# These functions allows us to concisely define the opt function below
	_tup = _x -> (optWrt == :traj ? (_x, params) : (traj, _x))
	x = (optWrt == :traj ? copy(traj) : copy(params))

	# Define the things needed for IPOPT
	g_L, g_U = gbounds(m, opt, traj, εs...)
	y0 = copy(traj[1:ny])
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, opt, _tup(x)..., y0)

	nnz = optWrt == :traj ? (opt.order == 1 ? Dgnnz(m, opt, traj) : DgnnzDirCol(m, opt, traj)) : N*ny*np
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = optWrt == :traj ? (opt.order == 1 ? Dgsparse!(rows, cols, values, m, opt, _tup(x)..., mode, ny, nu, N, δt) : DgsparseDirCol!(rows, cols, values, m, opt, _tup(x)..., mode, ny, nu, N, δt)) : Dgpsparse!(rows, cols, values, m, opt, _tup(x)..., mode, ny, nu, N, δt)
	x_L, x_U = optWrt == :traj ? xbounds(m, opt, N) : plimits(m)

	function eval_f(x::AbstractArray)
		_ro = robj(m, opt, _tup(x)...)
		return (_ro ⋅ _ro)
	end
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = ForwardDiff.gradient!(grad_f, eval_f, x)

	# # TEST the interface
	# println("TESTING")
	# println("xL", x_L, "xU", x_U)
	# println("gL", g_L, "gU", g_U)
	# g = similar(g_L)
	# eval_g(traj, g)
	# println("g", g)
	# println("f", eval_f(traj))
	# Df = similar(traj)
	# eval_grad_f(traj, Df)
	# println("Df", Df)

	# Create IPOPT problem
	prob = Ipopt.createProblem(
		length(x), # Number of variables
		x_L, # Variable lower bounds
		x_U, # Variable upper bounds
		length(g_L), # Number of constraints
		g_L,       # Constraint lower bounds
		g_U,       # Constraint upper bounds
		nnz,  # Number of non-zeros in Jacobian
		0,             # Number of non-zeros in Hessian
		eval_f,                     # Callback: objective function
		eval_g,                     # Callback: constraint evaluation
		eval_grad_f,                # Callback: objective function gradient
		eval_jac_g,                 # Callback: Jacobian evaluation
		nothing           # Callback: Hessian evaluation
	)
	Ipopt.addOption(prob, "hessian_approximation", "limited-memory")

	# Add options using kwargs
	for (k,v) in pairs(kwargs)
		println(k, " => ", v)
		Ipopt.addOption(prob, string(k), v)
	end

	# TODO: this should be an update only without need to setup. would need to update params.
	prob.x = x
	status = Ipopt.solveProblem(prob)

	if Ipopt.ApplicationReturnStatus[status] == :Infeasible_Problem_Detected
		# println("HIHI", prob.x)
	end

	# return status
	return prob
end

