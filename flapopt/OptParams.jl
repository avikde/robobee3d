
include("OptBase.jl") #< including this helps vscode reference the functions in there
using ProgressMeter # temp
#=========================================================================
Param opt
=========================================================================#

"Figure out the best δx direction to move the traj"
function paramδx(m::Model, opt::OptOptions, traj::AbstractArray, param0::AbstractArray, mult_x_L::AbstractArray, mult_x_U::AbstractArray)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# # Desired δx: lower u used
	# δx = copy(-0.1 * traj)# negative of the currently applied force
	# fill!(δx[1:(N+1)*ny], 0.0)

	# step in the direction of the active constraints
	δx = mult_x_L - mult_x_U

	return δx
end

"""
Set up a QP param step problem.
	δp is the decision var.
	P is set as a small regularization.
	The linear constraint is given by the IFT condition keeping δg = 0.

"""
function paramoptQPSetup(m::Model, opt::OptOptions, traj::AbstractArray; Preg=1e-3, settings...)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# Get number of constraints ng
	gL, gU = gbounds(m, opt, traj)
	ng = length(gL)
	np = pdims(m)
	
	# P and A are sparse matrices of type SparseMatrixCSC
	# TODO: A for param box constraints
	# # pick worst case for A = δg/δp. IC or symm constraints have no p; but assume Pk = δt δf/δp are all full
	# row = repeat(ny+1:(N+1)*ny, outer=np)
	# col = repeat(1:np, inner=N*ny)
	# val = ones(length(row))

	# everything will be updated - just need to set the sparsity
	mo = OSQP.Model()
	# OSQP.setup!(mo; P=sparse(ones(np,np)), q=ones(np), A=sparse(row, col, val, ng, np), l=ones(ng), u=ones(ng), settings...)
	OSQP.setup!(mo; P=sparse(ones(np,np)), q=ones(np), settings...) # no constraint for now
	return mo
end

"""
Q = prioritization matrix of size ? x nc, where nc is the size of g
"""
function paramopt(mo::Union{Nothing, OSQP.Model}, m::Model, opt::OptOptions, traj::AbstractArray, param0::AbstractArray, δx::AbstractArray, εs, Q::Union{Nothing, AbstractArray}; step=0.05, penalty=1e2)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	
	# NOTE: this is really more suited to the custom solver where Dg is already computed
	g_L, g_U = gbounds(m, opt, traj, εs...)
	Nc = length(g_L)
	Nx = length(traj)
	# Get Dg using the IPOPT functions
	nnz = Dgnnz(m, opt, traj)
	row = zeros(Int32, nnz)
	col = zeros(Int32, nnz)
	value = zeros(nnz)
	Dgsparse!(row, col, value, m, opt, traj, param0, :Structure, ny, nu, N, δt)
	Dgsparse!(row, col, value, m, opt, traj, param0, :Values, ny, nu, N, δt)
	# use sparse matrix from the COO format (same as used by IPOPT)
	Dg = sparse(row, col, value)

	# Gradient wrt params TODO: more insight than autograd?
	function gofp(pp)
		gg = Array{Any,1}(undef, length(g_L))
		gvalues!(gg, m, opt, traj, pp, traj[1:ny])
		return gg
	end
	dg_dp = convert(Array{Float64}, ForwardDiff.jacobian(gofp, param0))

	if isnothing(mo)
		# Non-QP version first
		if isnothing(Q)
			Q = I
		end
		δp = -(Q * dg_dp) \ (Q * Dg * δx)
		return param0 + step * δp
	else
		np = pdims(m)
		# QP version
		# l = -Dg * δx

		P = penalty * dg_dp' * dg_dp
		Px_new = zeros(np*(np+1)÷2)
		offs = 0
		for j = 1:np
			for i = 1:j
				offs += 1
				Px_new[offs] = P[i,j]
			end
		end

		# # Add a term related to the objective
		# function Jo(p)
		# 	_ro = robj(m, opt, traj, p)
		# 	return (_ro ⋅ _ro)
		# end
		# dJo_dp = convert(Array{Float64}, ForwardDiff.gradient(Jo, param0))
		
		# in the QP solution the "step size" only applies to the δx desired
		q = penalty * dg_dp' * Dg * δx * step #+ dJo_dp

		# println(Nn, size(Nn))
		OSQP.update!(mo, Px=Px_new; q=q)
		res = OSQP.solve!(mo)
		# print(fieldnames(typeof(res)))
		# error("HIHI")

		return param0 + res.x
	end
end

function paramoptJ(m::Model, opt::OptOptions, traj::AbstractArray, params0::AbstractArray, εs; step=0.05)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	
	function eval_f(x::AbstractArray)
		_ro = robj(m, opt, traj, x)
		return (_ro ⋅ _ro)
	end
	dJ_dp = ForwardDiff.gradient(eval_f, params0)

	return params0 - dJ_dp * step
end

# -----------------------------

function optboth(mo::Union{Nothing, OSQP.Model}, m::Model, opt::OptOptions, traj0::AbstractArray, param0::AbstractArray, εs, Q::Union{Nothing, AbstractMatrix}; step=0.05, penalty=1e2)
	prob = ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1, nlp_scaling_method="none")
	traj1 = prob.x
	# # with my modification to Ha/Coros g-preferred param opt
	# δx = paramδx(m, opt, traj1, param0, prob.mult_x_L, prob.mult_x_U)
	# param1 = paramopt(nothing, m, opt, traj1, param0, δx, εs, Q; step=step, penalty=penalty)

	# Test specific strategy
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj1)
	intu = sum([traj1[liu[1,k]] for k=1:N])
	intx = sum([traj1[liy[1,k]] for k=1:N])
	param1 = param0 + [-0.5 * intu / intx]

	return traj1, param1
end

# -------------

"For comparison; a dumb line search for the param"
function optnaive(mo::Union{Nothing, OSQP.Model}, m::Model, opt::OptOptions, traj0::AbstractArray, εs)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj0)

	ktest = collect(1:1.0:20)
	os = similar(ktest)

	@showprogress 1 "optnaive " for i = 1:length(ktest)
		param = [ktest[i]]
		prob = ipoptsolve(m, opt, traj0, param, εs, :traj; print_level=1, nlp_scaling_method="none")
		traj1 = prob.x
		os[i] = norm(traj1[(N+1)*ny+1:end]) # input force required
	end

	return ktest, os
end

# --------------
"Implement this"
function paramAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, R::AbstractArray)

end

function optAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, R::AbstractArray)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)

	paramAffine(m, opt, traj, param, R)
end

