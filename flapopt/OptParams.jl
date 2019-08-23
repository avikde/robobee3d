
include("OptBase.jl") #< including this helps vscode reference the functions in there

#=========================================================================
Param opt
=========================================================================#

"Figure out the best δx direction to move the traj"
function paramδx(m::Model, opt::OptOptions, traj::AbstractArray, param0::AbstractArray, mult_x_L::AbstractArray, mult_x_U::AbstractArray)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	# # Desired δx: lower u used
	# δx = copy(traj)
	# fill!(δx[1:(N+1)*ny], 0.0)
	# δx[(N+1)*ny+1:(N+1)*ny+N*nu] .= -δx[(N+1)*ny+1:(N+1)*ny+N*nu] # negative of the currently applied force

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
	# pick worst case for A = δg/δp. IC or symm constraints have no p; but assume Pk = δt δf/δp are all full
	row = repeat(ny+1:(N+1)*ny, outer=np)
	col = repeat(1:np, inner=N*ny)
	val = ones(length(row))

	# everything will be updated - just need to set the sparsity
	mo = OSQP.Model()
	OSQP.setup!(mo; P=spdiagm(0 => ones(np)), q=ones(np), A=sparse(row, col, val, ng, np), l=ones(ng), u=ones(ng), settings...)
	return mo
end

function paramopt(mo::Union{Nothing, OSQP.Model}, m::Model, opt::OptOptions, traj::AbstractArray, param0::AbstractArray, δx::AbstractArray, εs; step=0.05)
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
	dg_dp = ForwardDiff.jacobian(gofp, param0)

	if isnothing(mo)
		# Non-QP version first
		δp = -dg_dp \ Dg * δx
		return param0 + step * δp
	else
		np = pdims(m)
		# QP version
		q = zeros(np) # TODO:
		l = -Dg * δx

		Nn=nullspace(convert(Array{Float64}, dg_dp))
		println(Nn, size(Nn))
		# OSQP.update!(mo, Ax=; q=q, l=l, u=l)
		error("HIHI")
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
