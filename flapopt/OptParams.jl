
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
function paramAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, R::Tuple, scaleTraj=1.0)
	error("Implement this!")
end

"Implement this"
paramLumped(m::Model, param::AbstractArray) = error("Implement this")

"Return the transmission matrix that should multiply the state to go from OUTPUT to ACTUATOR"
function TmapAtoO(m::Model, T)
	error("Not implemented")
end

# lumped parameter vector
function getpt(m::Model, p)
	pb, Tarr = paramLumped(m, p)
	τ1, τ2 = Tarr
	# Nonlinear transmission: see https://github.com/avikde/robobee3d/pull/92
	return [pb*τ1; pb*τ2/τ1^2; 1/τ1; τ2/τ1^4], Tarr
end

"Override this"
function transmission(m::Model, y::AbstractArray, _param::Vector)
	yo = y
	T = 1.0
	τfun = x -> x
	τifun = x -> x
    return yo, T, τfun, τifun
end

"Helper function to reconstruct the traj (also test it)"
function reconstructTrajFromΔy(m::Model, opt::OptOptions, traj::AbstractArray, yo, Hk, B, Δy, pnew; test::Bool=false)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	np = length(pnew)
	Δyk = k -> Δy[(k-1)*ny+1 : k*ny]

	# Also convert the output traj with the Δy, new T, and inputs
	traj2 = copy(traj)
	# Calculate the new traj (which is in act coordinates, so needs scaling by T)
	ptnew, Tnew = getpt(m, pnew)
	for k=1:N+1
		# Go from output to act coords
		traj2[liy[:,k]] = (1.0 ./ TmapAtoO(m, Tnew)) .* (yo(k) + Δyk(k))
	end
	# Calculate the new inputs
	traj2[(N+1)*ny+1:end] = vcat([Tnew / δt * B' * Hk(k, Δyk(k), Δyk(k+1)) * ptnew for k=1:N]...) # compare to the "test" equation above

	# println("Test that the lower rows dyn constraint worked ", vcat([Bperp * Hk(k, Δyk(k), Δyk(k+1)) * ptnew for k=1:N]...) - eval_g_ret(prob.x))

	return traj2
end

"Mode=1 => opt, mode=2 ID. Fext(p) or hold constant"
function optAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, mode::Int, R::Tuple, εunact, plimsL, plimsU, scaleTraj=1.0; Fext_pdep::Bool=false, test=false, testTrajReconstruction=false, kwargs...)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
	np = length(param)
    ptTEST, TTEST = getpt(m, param) # NOTE the actual param values are only needed for the test mode
	npt = length(ptTEST)
    
	# Weights
	Ryy, Ryu, Ruu = R # NOTE Ryu is just weight on mech. power

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, R, scaleTraj; Fext_pdep=Fext_pdep)

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
	nx = np + (N+1)*ny # p,Δy

	# FIXME: print this out for now
	println("yomax = ", norm([yo(k)[1] for k=1:N], Inf))

	xlimsL = -1000 * ones(nx)
	xlimsU = 1000 * ones(nx)
	xlimsL[1:np] = plimsL
	xlimsU[1:np] = plimsU
	
	# ------------ Constraint: Bperp' * H(y + Δy) * pt is small enough (unactuated DOFs) -----------------
	nact = size(B, 2)
	nck = nq - nact # number of constraints for each k = # of unactuated DOFs ( = nunact)
	Bperp = (I - B*B')[nact+1:end,:] # s.t. Bperp*B = 0
	nc = N * nck# + np

	eval_g_pieces(k, Δyk, Δykp1, p) = Bperp * Hk(k, Δyk, Δykp1) * (getpt(m, p)[1])
	function eval_g_ret(x)
		Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
		# g .= x # OLD: 
		return vcat([eval_g_pieces(k, Δyk(k), Δyk(k+1), x[1:np]) for k=1:N]...)
	end
	# g1 = Array{Any,1}(undef, nc)
	eval_g(x::Vector, g::Vector) = g .= eval_g_ret(x)

	# ----------- Constraint Jac ----------------------------
	# Exploit sparsity in the nc*nx matrix. Each constraint depends on Δyk(k), Δyk(k+1), p
	Dgnnz = nc * (2*ny + np)

	# Function for IPOPT
	function eval_jac_g(x, imode, row::Vector{Int32}, col::Vector{Int32}, value)
		offs = 1
			
		if imode != :Structure
			Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
			p = x[1:np]

			for k=1:N
				# Now assemble the pieces
				dyk = ForwardDiff.jacobian(yy -> eval_g_pieces(k, yy, Δyk(k+1), p), Δyk(k))
				for i=1:nck
					for j=1:ny
						value[offs] = dyk[i,j]
						offs += 1
					end
				end
				dykp1 = ForwardDiff.jacobian(yy -> eval_g_pieces(k, Δyk(k), yy, p), Δyk(k+1))
				for i=1:nck
					for j=1:ny
						value[offs] = dykp1[i,j]
						offs += 1
					end
				end
				dp = ForwardDiff.jacobian(yy -> eval_g_pieces(k, Δyk(k), Δyk(k+1), yy), p)
				for i=1:nck
					for j=1:np
						value[offs] = dp[i,j]
						offs += 1
					end
				end
			end
		else
			for k=1:N
				for i=1:nck
					for j=1:ny
						row[offs] = (k-1)*nck + i
						col[offs] = np + (k-1)*ny + j
						offs += 1
					end
				end
				for i=1:nck
					for j=1:ny
						row[offs] = (k-1)*nck + i
						col[offs] = np + (k)*ny + j
						offs += 1
					end
				end
				for i=1:nck
					for j=1:np
						row[offs] = (k-1)*nck + i
						col[offs] = j
						offs += 1
					end
				end
			end
		end
	end

	glimsL = -εunact*ones(nc)
	glimsU = εunact*ones(nc)

	# ----------------------- Objective --------------------------------
    # If test is true, it will test the affine relation
    if test
        Hpb = zeros(nq, N)
		Bu = similar(Hpb)
		for k=1:N
			Hpb[:,k] = Hk(k, zeros(ny), zeros(ny)) * ptTEST
			Bu[:,k] = δt * B * umeas(k)[1]
		end
        display(Hpb - Bu)
        error("Tested")
	end
	
	# TODO: this is NOT INCLUDING the Δy in the calculation of u. Including these was resulting in a lot of IPOPT iterations and reconstruction failed -- need to investigate why. https://github.com/avikde/robobee3d/pull/80#issuecomment-541350179
	Quu = zeros(npt, npt)
	qyu = zeros(npt)
	qyy = 0
	
	# Matrices that contain sum of running actuator cost
	for k=1:N
		Hh = Hk(k, zeros(ny), zeros(ny)) #Hk(k, Δyk(k), Δyk(k+1))
		yok = yo(k)# + Δyk(k)
		if mode == 1
			Quu += Hh' * Ruu * Hh
			# Need output coords
			qyu += Ryu * (Hh' * [B * B'  zeros(ny-nq, ny-nq)] * yok)
			qyy += yok' * Ryy * yok # qyy * T^(-2)
		elseif mode == 2
			# Need to consider the unactuated rows too
			Quu += Hh' * Ruu * Hh
			# For ID, need uk
			qyu += (-Hh' * Ruu * (δt * B * umeas(k)))
		end
	end
	
	# See eval_f for how these are used to form the objective
	function eval_f(x)
		pt, Tarr = getpt(m, x[1:np])
		T = Tarr[1] # FIXME:
		# Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]

		# min Δy
		J = dot(x[np+1 : end], x[np+1 : end])

		if mode == 1
			# Total cost: 1/2 (T*pt)' * Quu * (T*pt) + 1/2 qyy * T^(-2) + qyu' * pt
			# TODO: quadratic version. for now just nonlinear
			J += 1/2 * ((T*pt)' * Quu * (T*pt) + qyy * T^(-2)) + qyu' * pt
		elseif mode == 2
			J += 1/2 * ((T*pt)' * Quu * (T*pt)) + qyu' * (T*pt)
		end

		return J
	end
	eval_grad_f(x, grad_f) = ForwardDiff.gradient!(grad_f, eval_f, x)

	# # Plot
	# display(Quu)
	# display(qyu)
	# Ts = collect(1.0:1.0:200.0)
	# fs = [eval_f([param[1], T]) for T in Ts]
    # plot(Ts, fs, xlabel="T")
    # gui()
	
	# Create IPOPT problem
	prob = Ipopt.createProblem(
		nx, # Number of variables
		xlimsL, # Variable lower bounds
		xlimsU, # Variable upper bounds
		nc, # Number of constraints
		glimsL,       # Constraint lower bounds
		glimsU,       # Constraint upper bounds
		Dgnnz,  # Number of non-zeros in Jacobian
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

	# Solve
	prob.x = [copy(param); zeros(nx - np)]
	status = Ipopt.solveProblem(prob)
	pnew = prob.x[1:np]
	trajnew = reconstructTrajFromΔy(m, opt, traj, yo, Hk, B, prob.x[np+1:end], pnew)
	unactErr = eval_g_ret(prob.x)

	if testTrajReconstruction
		# Test traj reconstruction:
		Hk, yo, umeas, B, N = paramAffine(m, opt, trajnew, pnew, R; Fext_pdep=Fext_pdep)
		eval_g_ret2(p) = vcat([Bperp * Hk(k, zeros(ny), zeros(ny)) * (getpt(m, p)[1]) for k=1:N]...)
		display(unactErr')
		display(eval_g_ret2(pnew)')
		error("Tested")
	end

	return pnew, eval_f, trajnew, unactErr

	# # Without that T, can just use OSQP -------------------------
	# mo = OSQP.Model()
	# # OSQP.setup!(mo; P=sparse(ones(np,np)), q=ones(np), A=sparse(row, col, val, ng, np), l=ones(ng), u=ones(ng), settings...)
	# OSQP.setup!(mo; P=sparse(Rp), q=zeros(np)) # no constraint for now

	# # OSQP.update!(mo; q=q)
	# res = OSQP.solve!(mo)

	# return res.x
end

