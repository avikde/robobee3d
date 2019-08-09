
include("Model.jl") #< including this helps vscode reference the functions in there

struct OptWorkspace
	x::Vector
	g::Vector
	# TODO: sparse matrices for these spzeros
	Dg::Matrix
	DgTg::Vector
	∇J::Vector
	HJ::Matrix
	OptWorkspace(Nx::Int, Nc::Int) = new(
		zeros(Nx), 
		zeros(Nc), 
		zeros(Nc, Nx), 
		zeros(Nx), 
		zeros(Nx), 
		zeros(Nx, Nx)
	)
end

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
function gvalues!(gout::Vector{T}, m::Model, opt::OptOptions, traj::Vector{T}, params::Vector{T}, y0::AbstractArray{T}) where {T}
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	li = LinearIndices((1:ny, 1:(N+2))) # to N+2

	yk = k -> @view traj[liy[:,k]]
	uk = k -> @view traj[liu[:,k]]

	gout[1:(N+1)*ny] = -(@view traj[1:(N+1)*ny])
	# Dynamics constraint
	for k = 1:N
		gout[li[:,k+1]] += ddynamics(m, yk(k), uk(k), params, δt)
	end

	# Initial condition
	gout[li[:,1]] = y0 - yk(1)

	# Periodicity or symmetry
	if opt.boundaryConstraint == SYMMETRIC
		# FIXME: for now hardcoded symmetry G(y) = -y
		gout[li[:,N+2]] = -yk(1) - yk(N+1)
	end

	return
end
#=
# "Bounds corresponding to the constraints above"
function gbounds(m::Model, opt::OptOptions, traj::Vector)::Tuple{Vector, Vector}
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)

	# println("CALLED gbounds with $(ny) $(nu) $(N) $(pointer_from_objref(traj))")
	g_L = [-traj[1:ny]; zeros(N*ny)]
    g_U = [-traj[1:ny]; zeros(N*ny)]
    # first N*ny = 0 (dynamics)
    return g_L, g_U
end

function Dgnnz(m::Model, opt::OptOptions, N::Int)::Int
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# Assuming the Jacobians are dense. The terms below correspond to the "-I"s, the "I + δt A"'s, the "δt B"'s
	return ny*(N+1) + ny^2*N + ny*nu*N
end

function Dgsparse!(row::Vector{Int32}, col::Vector{Int32}, value::Vector, m::Model, opt::OptOptions, traj::Vector, params::Vector, mode)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)

	# NOTE: traj is NULL when in :Structure mode
	if mode != :Structure
		δt = opt.vart ? traj[end] : opt.fixedδt
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
			dlin!(Ak, Bk, m, traj[@view liy[:,k]], traj[@view liu[:,k]], params, δt)
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
	return
end
=#
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

@enum OptVar WRT_TRAJ WRT_PARAMS

"""Custom solver"""
function csSolve!(wk::OptWorkspace, m::Model, opt::OptOptions, traj0::AbstractArray{T}, params0::AbstractArray, optWrt::OptVar; μs::Array{Float64}=[1e-1], Ninner::Int=1) where {T}
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj0)

	# These functions allows us to concisely define the opt function below
	_tup = _x -> (optWrt == WRT_TRAJ ? (_x, params0) : (traj0, _x))
	robjx = _x -> robj(m, opt, _tup(_x)...)
	x = (optWrt == WRT_TRAJ ? copy(traj0) : copy(params0))

	trajp = _tup(x)
	# functions for views
	yk = k -> @view trajp[1][liy[:,k]]
	uk = k -> @view trajp[1][liu[:,k]]
	gk = k -> @view wk.g[liy[:,k]]

	# get (y,u,p,δt) at time k--point at which to evaluate dynamics
	_yupt = k::Int -> (yk(k), uk(k), trajp[2], δt)

	# Constraint bounds
	x_L, x_U = optWrt == WRT_TRAJ ? xbounds(m, opt, N) : (fill(-Inf, size(params0)), fill(Inf, size(params0)))

	if optWrt == WRT_TRAJ
		for i = 1:ny*(N+1)
			wk.Dg[i,i] = -1.0
		end
		if opt.boundaryConstraint == SYMMETRIC
			for j = 1:ny
				wk.Dg[(N+1) * ny + j, j] = -1.0
				wk.Dg[(N+1) * ny + j, (N) * ny + j] = -1.0
			end
		end
		Ak = zeros(ny, ny)
		Bk = zeros(ny, nu)
	elseif optWrt == WRT_PARAMS
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
				return (_ro ⋅ _ro) + μ/2 * ((wk.g ⋅ wk.g) + sum(Ψ.(_x - x_U) + Ψ.(x_L - _x)))
			end

			# Compute Jacobian and Hessian
			for k = 1:N
				if optWrt == WRT_TRAJ
					dlin!(Ak, Bk, m, _yupt(k)...)
					# Dg_y' * g = [-g0 + A1^T g1, ..., -g(N-1) + AN^T gN, -gN]
					wk.DgTg[liy[:,k]] = -gk(k) + Ak' * gk(k+1)
					# Dg_u' * g = [B1^T g1, ..., BN^T gN]
					wk.DgTg[liu[:,k]] = Bk' * gk(k+1)
					# Dg_δt' * g = 0

					# TODO: better Dg' Dg computation that doesn't compute Dg
					wk.Dg[liy[:,k+1], liy[:,k]] = Ak
					wk.Dg[liy[:,k+1], liu[:,k]] = Bk
				elseif optWrt == WRT_PARAMS
					dlinp!(Pk, m, _yupt(k)...)
					# Dg0 = 0
					wk.Dg[liy[:,k+1], :] = Pk
					wk.DgTg = wk.DgTg + Pk' * gk(k)
				end
			end
			
			if optWrt == WRT_TRAJ
				wk.DgTg[liy[:,N+1]] = -gk(N+1)
			end

			# No special structure for the objective, but we only need the gradient and no Hessian
			ro = robjx(x)
			Dro = ForwardDiff.jacobian(robjx, x)

			# Gradient: most terms are quadratic and hence use the Gauss-Newton approx; otherwise use the ineq constraint and its special "diagonal" form
			x_Udiff = x - x_U
			x_Ldiff = x_L - x
			wk.∇J[:] = Dro' * ro + μ * (wk.DgTg + 1/2 * (dΨ.(x_Udiff) - dΨ.(x_Ldiff)))
			wk.HJ[:] = Dro' * Dro + μ * (wk.Dg' * wk.Dg)
			wk.HJ[diagind(wk.HJ)] += μ/2 * (ddΨ.(x_Udiff) + ddΨ.(x_Ldiff))
			# Regularization, and then we know this is pos def
			wk.HJ[diagind(wk.HJ)] .+= opt.hessReg

			# # Gradient descent
			# v .= -∇J
			# Newton or Gauss-Newton. Use PositiveFactorizations.jl to ensure psd Hessian
			v .= -wk.HJ\wk.∇J #-(cholesky(Positive, wk.HJ) \ wk.∇J)

			J0 = Jx(x)
			J1 = csBacktrackingLineSearch!(x1, x, wk.∇J, v, J0, Jx; α=0.2, β=0.7)
			x .= x1
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

function csAlternateSolve(m::Model, opt::OptOptions, traj0::Vector, params0::Vector, NaltSteps::Int=1; μst::Array{Float64}=[1e-1], Ninnert::Int=1, μsp::Array{Float64}=[1e-1], Ninnerp::Int=1)
	# reshape into Nx1 matrices
	trajs = reshape(copy(traj0), :, 1)
	params = reshape(copy(params0), :, 1)
	# Append columns for each step
	for isteps = 1:NaltSteps
		@time trajs = [trajs csSolve(m, opt, trajs[:,end], params[:,end], WRT_TRAJ; Ninner=Ninnert, μs=μst)]
		@time params = [params csSolve(m, opt, trajs[:,end], params[:,end], WRT_PARAMS; Ninner=Ninnerp, μs=μsp)]
	end
	return trajs, params
end

#=========================================================================
Solver interface
=========================================================================#

function nloptsetup(m::Model, opt::OptOptions, traj::Vector, params::Vector)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)

	# Define the things needed for IPOPT
	x_L, x_U = xbounds(m, opt, N)
	g_L, g_U = gbounds(m, opt, traj)
	eval_g(x::Vector, g::Vector) = gvalues!(g, m, opt, x, params, traj[1:ny])
	eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector) = Dgsparse!(rows, cols, values, m, opt, x, params, mode)
	eval_f(x::Vector{Float64}) = Jobj(m, opt, x, params)
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = ForwardDiff.gradient!(grad_f, eval_f, x)

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
