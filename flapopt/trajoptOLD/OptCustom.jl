
include("OptBase.jl") #< including this helps vscode reference the functions in there

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
	np = pdims(m)

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
		Pk = zeros(ny, np)
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
					wk.DgTg .= wk.DgTg + Pk' * gk(k)
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

			# Gradient descent
			# v .= -wk.∇J
			# Newton or Gauss-Newton. Use PositiveFactorizations.jl to ensure psd Hessian
			# v .= -(cholesky(Positive, wk.HJ) \ wk.∇J)# #
			v .= -wk.HJ\wk.∇J

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

csAlwaysFeasible(x::Vector) = true

function csBacktrackingLineSearch!(x1::Vector, x0::Vector, ∇J0::Vector, v::Vector, J0::Float64, Jcallable; α::Float64=0.45, β::Float64=0.9, isFeasible::Function=csAlwaysFeasible)
	σ = 1
	# search for step size
	while true
		σ = β * σ
		x1 .= x0 + σ * v
		J1 = Jcallable(x1)
		# debug line search
		# println("J0=$(round(J0; sigdigits=4)), J1=$(round(J1; sigdigits=4)), σ=$(round(σ; sigdigits=6))")
		if (J1 < J0 + α * σ * ∇J0' * v && isFeasible(x1)) || σ < 1e-6
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
	wkt = OptWorkspace(Ntraj(m, opt, N), (N+2)*ny)
	# np = length(params0)
	# wkp = OptWorkspace(np, (N+2)*ny)

	# Append columns for each step
	for isteps = 1:NaltSteps
		@time trajs = [trajs csSolve!(wkt, m, opt, trajs[:,end], params[:,end], :traj; Ninner=Ninnert, μs=μst)]
		# @time params = [params csSolve!(wkp, m, opt, trajs[:,end], params[:,end], WRT_PARAMS; Ninner=Ninnerp, μs=μsp)]
	end
	
	return trajs, params, wkt
end
