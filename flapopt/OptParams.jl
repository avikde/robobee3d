
include("OptBase.jl") #< including this helps vscode reference the functions in there
include("OptCustom.jl")
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
function paramAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, R::Tuple)
	error("Implement this!")
end

"Implement this"
paramLumped(m::Model, param::AbstractArray) = error("Implement this")

"Mode=1 => opt, mode=2 ID"
function optAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, mode::Int, R::Tuple; hessreg::Float64=0, kwargs...)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
    nq = ny÷2
    # lumped parameter vector
    pb, T = paramLumped(m, param) # NOTE the actual param values are only needed for the test mode
    # pb = [T^2, cbar*T]
    pt = [pb; T^(-2)]
    npt = length(pt) # add T^(-2)
    
    # Weights
    Ryy, Ryu, Ruu = R # NOTE Ryu is just weight on mech. power

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, R)

    # If test is true, it will test the affine relation
    test = true
    if test
        Hpb = zeros(nq, N)
        Bu = similar(Hpb)
    end
    
    # Quu = zeros(npt, npt)
    # qyu = zeros(npt)
    # qyumeas = zeros(npt) # for ID
    # qyy = 0
    for k=1:N
        # Quu += Hh' * B * Ruu * B' * Hh # (T*pt)' * Quu * (T*pt)
        # # Need output coords
        # yo = [T, 1.0, T, 1.0] .* yk(k)
        # qyu += Ryu * (Hh' * [B * B'  zeros(2, 2)] * yo) # qyu' * pt
        # qyy += yo' * Ryy * yo # qyy * T^(-2)
        # # For ID, need uk
        # qyumeas -= Hh' * B * Ruu * (δt * uk(k))

        if test
            Hpb[:,k] = Hk(k) * pt
            Bu[:,k] = δt * B * umeas(k)[1] / T
        end
    end
    if test
        display(Hpb - Bu)
        error("Tested")
    # else
    #     return Quu, qyu, qyy, qyumeas
    end

	error("HERE")

	# GN -------------------------
	# function pFeasible(p)
	# 	cbar, T = p
	# 	return T >= 1 # FIXME: testing
	# end
	# np = 2

	# # variable to update
	# x = copy(param)
	# x1 = similar(x)

	# # Gauss-Newton iterations with the feasible space of p
	# for stepi=1:10
	# 	# Cost function for this step
	# 	Jx = p -> pb(p)' * Rp * pb(p)
	# 	rx = p -> S * pb(p)
	# 	# Current jacobian
	# 	r0 = rx(x)
	# 	Dr0 = ForwardDiff.jacobian(rx, x)
	# 	# Gauss-Newton
	# 	∇J = Dr0' * r0
	# 	HJ = Dr0' * Dr0 + hessreg * I
	# 	v = -HJ\∇J #descent direction

	# 	J0 = Jx(x)
	# 	J1 = csBacktrackingLineSearch!(x1, x, ∇J, v, J0, Jx; α=0.2, β=0.7, isFeasible=pFeasible)
	# 	x .= x1
	# end
	
	# return x

	# IPOPT ---------------------------
	
	Aconstraint = [1.0  1.0] # testing linear constraint - useless for now
	eval_g(x::Vector, g::Vector) = (g .= Aconstraint * x)
	function eval_jac_g(x::Vector{Float64}, mode, row::Vector{Int32}, col::Vector{Int32}, value::Vector)
		if mode != :Structure
			value[1] = Aconstraint[1]
			value[2] = Aconstraint[2]
		else
			row[1] = row[2] = 1
			col[1] = 1
			col[2] = 2
		end
	end

	function eval_f(x::AbstractArray)
		pb, T = paramLumped(m, x)
		pt = [pb; T^(-2)]
		if mode == 1
			# Total cost: 1/2 (T*pt)' * Quu * (T*pt) + 1/2 qyy * T^(-2) + qyu' * pt
			# TODO: quadratic version. for now just nonlinear
			return 1/2 * ((T*pt)' * Quu * (T*pt) + qyy * T^(-2)) + qyu' * pt
		elseif mode == 2
			return 1/2 * ((T*pt)' * Quu * (T*pt)) + qyumeas' * (T*pt)
		end
	end
	eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64}) = ForwardDiff.gradient!(grad_f, eval_f, x)

	# # Plot
	# display(Quu)
	# display(qyu)
	# Ts = collect(1.0:1.0:200.0)
	# fs = [eval_f([param[1], T]) for T in Ts]
    # plot(Ts, fs, xlabel="T")
    # gui()
	
	# Create IPOPT problem
	prob = Ipopt.createProblem(
		length(param), # Number of variables
		[0.1, 0.1], # Variable lower bounds
		[10.0, 1000.0], # Variable upper bounds
		1, # Number of constraints
		[-1000.0],       # Constraint lower bounds
		[1000.0],       # Constraint upper bounds
		2,  # Number of non-zeros in Jacobian
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
		# println(k, " => ", v)
		Ipopt.addOption(prob, string(k), v)
	end

	# TODO: this should be an update only without need to setup. would need to update params.
	prob.x = copy(param)
	status = Ipopt.solveProblem(prob)
	return prob.x, eval_f

	# # Without that T, can just use OSQP -------------------------
	# mo = OSQP.Model()
	# # OSQP.setup!(mo; P=sparse(ones(np,np)), q=ones(np), A=sparse(row, col, val, ng, np), l=ones(ng), u=ones(ng), settings...)
	# OSQP.setup!(mo; P=sparse(Rp), q=zeros(np)) # no constraint for now

	# # OSQP.update!(mo; q=q)
	# res = OSQP.solve!(mo)

	# return res.x
end

