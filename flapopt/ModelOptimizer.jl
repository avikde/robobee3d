
include("Model.jl") #< including this helps vscode reference the functions in there

struct OptWorkspace
	x::Vector
	g::Vector
	# TODO: sparse matrices for these spzeros
	Dg::Matrix
	DgTg::Vector
	∇J::Vector
	HJ::Matrix
	λ::Vector # Augmented Lagrangian
	OptWorkspace(Nx::Int, Nc::Int) = new(
		zeros(Nx), 
		zeros(Nc), 
		zeros(Nc, Nx), 
		zeros(Nx), 
		zeros(Nx), 
		zeros(Nx, Nx), 
		zeros(Nc)
	)
end

#=========================================================================
Dynamics constraint
=========================================================================#

function gkDirCol(m::Model, params, yk, ykp1, uk, ukp1, δt, fk=nothing)
	uc = (uk + ukp1) / 2
	if isnothing(fk)
		fk = dydt(m, yk, uk, params)
	end
	fkp1 = dydt(m, ykp1, ukp1, params)
	# interpolate http://underactuated.mit.edu/underactuated.html?chapter=trajopt
	yc = 1/2 * (yk + ykp1) + δt/8 * (fk - fkp1)
	ycdotδt = -3/(2) * (yk - ykp1) - δt/4 * (fk + fkp1)
	return -ycdotδt + δt * dydt(m, yc, uc, params), fkp1
end

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
function gvalues!(gout::Vector{T}, m::Model, opt::OptOptions, traj::Vector{T}, params::Vector{T}, y0::AbstractArray{T}) where {T}
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	li = LinearIndices((1:ny, 1:(N+2))) # to N+2

	yk = k -> @view traj[liy[:,k]]
	uk = k -> @view traj[liu[:,k]]

	# Initial condition
	gout[li[:,1]] = y0 - yk(1)

	# Dynamics constraint
	if opt.order > 1
		fk = similar(y0)
		fkp1 = dydt(m, yk(1), uk(1), params) # initialize for usage below
		uc = zeros(nu)
		yc = similar(y0)
		ycdotδt = similar(y0)
		ukp1 = similar(uc)
	end

	for k = 1:N
		if opt.order == 1
			gout[li[:,k+1]] = -yk(k+1) + ddynamics(m, yk(k), uk(k), params, δt)
		else # assume collocation
			# Get collocation point states/inputs
			if k < N
				ukp1 .= uk(k+1)
			elseif k == N && opt.boundaryConstraint == :symmetric
				ukp1 .= -uk(1)
			else
				ukp1 .= uk(k) # don't have any more u's in the current parameterization
			end
			# Can reuse fkp1 as fk for the next loop iterate
			fk .= fkp1
			# collocation constraint
			gout[li[:,k+1]], fkp1 = gkDirCol(m, yk(k), yk(k+1), uk(k), ukp1, δt, fk)
		end
	end

	# Periodicity or symmetry
	if opt.boundaryConstraint == :symmetric
		# FIXME: for now hardcoded symmetry G(y) = -y
		gout[li[:,N+2]] = -yk(1) - yk(N+1)
	end

	return
end

#=========================================================================
Custom solver
=========================================================================#


"""Indicator function to use inequality constraints in a penalty method.
Following Geilinger et al (2018) skaterbots. Ref Bern et al (2017).

A scalar inequality f(x) <= b should be modeled as Ψ(f(x) - b) and added to the cost."""
function Ψ(x; ε::Float64=0.1)
    if x <= -ε
        return 0
	elseif x > -ε && x < ε
        return x^3/(6ε) + x^2/2 + ε*x/2 + ε^2/6
    else
		return x^2 + ε^2/3
	end
end

function dΨ(x; ε::Float64=0.1)
    if x <= -ε
        return 0
	elseif x > -ε && x < ε
        return x^2/(2ε) + x + ε/2
    else
		return 2x
	end
end

function ddΨ(x; ε::Float64=0.1)
    if x <= -ε
        return 0
	elseif x > -ε && x < ε
        return x/(ε) + 1
    else
		return 2
	end
end

"""Custom solver"""
function csSolve!(wk::OptWorkspace, m::Model, opt::OptOptions, traj0::AbstractArray, params0::AbstractArray, optWrt::Symbol; μs::Array{Float64}=[1e-1], Ninner::Int=1)
	optWrt in OptVar || throw(ArgumentError("invalid optWrt: $optWrt"))

	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj0)

	# These functions allows us to concisely define the opt function below
	_tup = _x -> (optWrt == :traj ? (_x, params0) : (traj0, _x))
	robjx = _x -> robj(m, opt, _tup(_x)...)
	x = (optWrt == :traj ? copy(traj0) : copy(params0))

	trajp = _tup(x)
	# functions for views
	yk = k -> @view trajp[1][liy[:,k]]
	uk = k -> @view trajp[1][liu[:,k]]
	gk = k -> @view wk.g[liy[:,k]]

	# get (y,u,p,δt) at time k--point at which to evaluate dynamics
	_yupt = k::Int -> (yk(k), uk(k), trajp[2], δt)

	# Constraint bounds
	x_L, x_U = optWrt == :traj ? xbounds(m, opt, N) : (fill(-Inf, size(params0)), fill(Inf, size(params0)))
	fill!(wk.λ, 0.0)

	if optWrt == :traj
		for i = 1:ny*(N+1)
			wk.Dg[i,i] = -1.0
		end
		if opt.boundaryConstraint == :symmetric
			for j = 1:ny
				wk.Dg[(N+1) * ny + j, j] = -1.0
				wk.Dg[(N+1) * ny + j, (N) * ny + j] = -1.0
			end
		end
		Ak = zeros(ny, ny)
		Bk = zeros(ny, nu)
	else
		Pk = zeros(ny, Nx)
	end
	
	# Take some number of steps
	# Preallocate output
	x1 = similar(x)
	v = similar(wk.∇J)

	for μ in μs
		for stepi = 1:Ninner
			# One step
			gvalues!(wk.g, m, opt, _tup(x)..., @view traj0[1:ny])

			# Cost function for this step
			function Jx(_x::AbstractArray)::Real
				gvalues!(wk.g, m, opt, _tup(_x)..., @view traj0[1:ny])
				_ro = robj(m, opt, _tup(_x)...)
				return (_ro ⋅ _ro) + μ/2 * ((wk.g ⋅ wk.g) + sum(Ψ.(_x - x_U) + Ψ.(x_L - _x))) + (wk.g ⋅ wk.λ)
			end

			# Compute Jacobian and Hessian
			for k = 1:N
				if optWrt == :traj
					dlin!(Ak, Bk, m, _yupt(k)...)
					# Dg_y' * g = [-g0 + A1^T g1, ..., -g(N-1) + AN^T gN, -gN]
					wk.DgTg[liy[:,k]] = -gk(k) + Ak' * gk(k+1)
					# Dg_u' * g = [B1^T g1, ..., BN^T gN]
					wk.DgTg[liu[:,k]] = Bk' * gk(k+1)
					# Dg_δt' * g = 0

					# TODO: better Dg' Dg computation that doesn't compute Dg
					wk.Dg[liy[:,k+1], liy[:,k]] = Ak
					wk.Dg[liy[:,k+1], liu[:,k]] = Bk
				else
					dlinp!(Pk, m, _yupt(k)...)
					# Dg0 = 0
					wk.Dg[liy[:,k+1], :] = Pk
					wk.DgTg = wk.DgTg + Pk' * gk(k)
				end
			end
			
			if optWrt == :traj
				wk.DgTg[liy[:,N+1]] = -gk(N+1)
			end

			# No special structure for the objective, but we only need the gradient and no Hessian
			ro = robjx(x)
			Dro = ForwardDiff.jacobian(robjx, x)

			# Gradient: most terms are quadratic and hence use the Gauss-Newton approx; otherwise use the ineq constraint and its special "diagonal" form
			x_Udiff = x - x_U
			x_Ldiff = x_L - x
			wk.∇J[:] = Dro' * ro + μ * (wk.DgTg + 1/2 * (dΨ.(x_Udiff) - dΨ.(x_Ldiff))) + (wk.Dg' * wk.λ)
			wk.HJ[:] = Dro' * Dro + μ * (wk.Dg' * wk.Dg)
			wk.HJ[diagind(wk.HJ)] += μ/2 * (ddΨ.(x_Udiff) + ddΨ.(x_Ldiff))
			# Regularization, and then we know this is pos def
			wk.HJ[diagind(wk.HJ)] .+= opt.hessReg

			# # Gradient descent
			# v .= -∇J
			# Newton or Gauss-Newton. Use PositiveFactorizations.jl to ensure psd Hessian
			# v = -(cholesky(Positive, wk.HJ) \ wk.∇J)# #
			v = -wk.HJ\wk.∇J

			J0 = Jx(x)
			J1 = csBacktrackingLineSearch!(x1, x, wk.∇J, v, J0, Jx; α=0.2, β=0.7)
			x .= x1
			# Update augmented Lagrangian
			gvalues!(wk.g, m, opt, _tup(x)..., @view traj0[1:ny])
			if opt.augLag
				wk.λ[:] -= μ/2 * wk.g
			end
			println("μ=$(μ)\tstep=$(stepi)\tJ $(round(J0;sigdigits=4)) → $(round(J1;sigdigits=4))")
		end
	end
	return x
end

function csBacktrackingLineSearch!(x1::Vector, x0::Vector, ∇J0::Vector, v::Vector, J0::Float64, Jcallable; α::Float64=0.45, β::Float64=0.9)
	σ = 1
	# search for step size
	while true
		σ = β * σ
		x1 .= x0 + σ * v
		J1 = Jcallable(x1)
		# debug line search
		# println("J0=$(round(J0; sigdigits=4)), J1=$(round(J1; sigdigits=4)), σ=$(round(σ; sigdigits=6))")
		if J1 < J0 + α * σ * ∇J0' * v || σ < 1e-6
			return J1
		end
	end
	return J0
end

function csAlternateSolve(m::Model, opt::OptOptions, traj0::AbstractArray, params0::AbstractArray, NaltSteps::Int=1; μst::Array{Float64}=[1e-1], Ninnert::Int=1, μsp::Array{Float64}=[1e-1], Ninnerp::Int=1)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj0)

	# reshape into Nx1 matrices
	trajs = reshape(copy(traj0), :, 1)
	params = reshape(copy(params0), :, 1)
	# Create workspaces
	wkt = OptWorkspace((N+1)*ny + N*nu + 1, (N+2)*ny)
	# np = length(params0)
	# wkp = OptWorkspace(np, (N+2)*ny)

	# Append columns for each step
	for isteps = 1:NaltSteps
		@time trajs = [trajs csSolve!(wkt, m, opt, trajs[:,end], params[:,end], :traj; Ninner=Ninnert, μs=μst)]
		# @time params = [params csSolve!(wkp, m, opt, trajs[:,end], params[:,end], WRT_PARAMS; Ninner=Ninnerp, μs=μsp)]
	end
	
	return trajs, params, wkt
end

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


function nloptsetup(m::Model, opt::OptOptions, traj::Vector, params::Vector; kwargs...)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	δt = opt.vart ? traj[end] : opt.fixedδt

	# Define the things needed for IPOPT
	x_L, x_U = xbounds(m, opt, N)
	g_L, g_U = gbounds(m, opt, traj, 100., 0.01, 0.)
	y0 = copy(traj[1:ny])
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, opt, x, params, y0)
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = Dgsparse!(rows, cols, values, m, opt, x, params, mode, ny, nu, N, δt)
	function eval_f(x::AbstractArray)
		_ro = robj(m, opt, x, params)
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
		length(traj), # Number of variables
		x_L, # Variable lower bounds
		x_U, # Variable upper bounds
		length(g_L), # Number of constraints
		g_L,       # Constraint lower bounds
		g_U,       # Constraint upper bounds
		Dgnnz(m, opt, traj),  # Number of non-zeros in Jacobian
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
	prob.x = copy(traj)

	return prob
end

function nloptsolve(prob)
	status = Ipopt.solveProblem(prob)

	if Ipopt.ApplicationReturnStatus[status] == :Infeasible_Problem_Detected
		# println("HIHI", prob.x)
	end

	return status
end

#=========================================================================
Visualization
=========================================================================#

# TODO: args...
function visualizeConstraintViolations(m::Model, opt::OptOptions, params::Vector, traj0::Vector, traj1::Vector)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj0)
	g_L, g_U = gbounds(m, opt, traj0)

	g0 = similar(g_L)
	g1 = similar(g0)
	gvalues!(g0, m, opt, traj0, params, traj0[1:ny])
	gvalues!(g1, m, opt, traj1, params, traj0[1:ny])

	println(g0, g1)
	pl2 = plot([g0,g1], marker=:auto, title="Constraint violations")
	hline!(pl2, [0], color="black", alpha=0.3)
end

