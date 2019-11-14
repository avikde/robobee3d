
include("OptBase.jl") #< including this helps vscode reference the functions in there

# --------------
"Implement this"
function paramAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, scaleTraj=1.0; Fext_pdep::Bool=true)
	error("Implement this!")
end

"Implement this"
paramLumped(m::Model, param::AbstractArray) = error("Implement this")

# lumped parameter vector
function getpt(m::Model, p)
	pb, Tarr = paramLumped(m, p)
	τ1, τ2 = Tarr
	# Nonlinear transmission: see https://github.com/avikde/robobee3d/pull/92
	return [pb*τ1; pb*τ2/τ1^2; 1/τ1; τ2/τ1^4], Tarr
end

"Override this. Input y can be in actuator or output coordinates. If o2a is false, maps A to O, else maps O to A"
function transmission(m::Model, y::AbstractArray, _param::Vector; o2a=false)
	y2 = y
	T = 1.0
	τfun = x -> x
	τifun = x -> x
    return y2, T, τfun, τifun
end

"σomax is an output strain limit. This is the only transmission constraint for now, but others can be added."
function gtransmission(m::Model, param, σomax)
	# Get both transmission coeffs
	pbb, Tarrr = paramLumped(m, param)
	τ1, τ2 = Tarrr
	return σomax/τ1 - σomax^3/3 * τ2/τ1^4
end

"Helper function to reconstruct the traj (also test it). trajAct true=>traj is in act coords (else output)"
function reconstructTrajFromΔy(m::Model, opt::OptOptions, traj::AbstractArray, yo, Hk, B, Δy, pnew, trajAct=true; test::Bool=false)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	np = length(pnew)
	Δyk = k -> Δy[(k-1)*ny+1 : k*ny]

	# Also convert the output traj with the Δy, new T, and inputs
	traj2 = copy(traj)
	# Calculate the new traj (which is in act coordinates, so needs scaling by T)
	ptnew, Tnew = getpt(m, pnew)
	for k=1:N+1
		if trajAct
			# Go from output to act coords
			ya, Tk = transmission(m, yo(k) + Δyk(k), pnew; o2a=true)[1:2]
			traj2[liy[:,k]] = ya
		else
			traj2[liy[:,k]] = yo(k) + Δyk(k)
		end
		# Calculate the new inputs
		if k <= N
			traj2[liu[:,k]] = 1 / δt * B' * Hk(k, Δyk(k), Δyk(k+1)) * ptnew # compare to the "test" equation above
		end
	end

	# println("Test that the lower rows dyn constraint worked ", vcat([Bperp * Hk(k, Δyk(k), Δyk(k+1)) * ptnew for k=1:N]...) - eval_g_ret(prob.x))

	return traj2
end

function affineTest(m, opt, traj, param)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
    ptTEST, TTEST = getpt(m, param) # NOTE the actual param values are only needed for the test mode

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, scaleTraj=1.0; Fext_pdep=Fext_pdep)
	Hpb = zeros(nq, N)
	Bu = similar(Hpb)
	for k=1:N
		Hpb[:,k] = Hk(k, zeros(ny), zeros(ny)) * ptTEST
		Bu[:,k] = δt * B * umeas(k)[1]
	end
	display(Hpb - Bu)
	error("Tested")
end

"Helper function for optAffine. See optAffine for the def of x"
function paramOptObjQuadratic(m, mode, np, npt, ny, nq, δt, Hk, yo, umeas, B, N, R)
	Ryy, Ryu, Ruu = R # NOTE Ryu is just weight on mech. power
	
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
	
	function eval_f(x)
		pt, Tarr = getpt(m, x[1:np])

		# min Δy
		J = dot(x[np+1 : end], x[np+1 : end])

		if mode == 1
			# FIXME: took out qyy component for nonlinear transmission
			J += 1/2 * (pt' * Quu * pt) + qyu' * pt# + qyy * T^(-2)) 
		elseif mode == 2
			J += 1/2 * (pt' * Quu * pt) + qyu' * pt
		end

		return J
	end
	eval_grad_f(x, grad_f) = ForwardDiff.gradient!(grad_f, eval_f, x)

	return eval_f, eval_grad_f
end

"Constraints for IPOPT: g = [gunact; gpolycon; gtransmission]
- gunact = unactuated error = #unactuated DOFS * N
- gpolycon = Cp * p <= dp. But the user can pass in (and default is) Cp 0xX => has no effect
- gtransmission: actuator strain limit.
"
function paramOptConstraint(m, mode, np, ny, nq, δt, Hk, yo, umeas, B, N, Cp, dp, σamax, σomax, εunact, τinds, uinfnorm; nonlinTransmission=true)
	# Unactuated constraint: Bperp' * H(y + Δy) * pt is small enough (unactuated DOFs) 
	nact = size(B, 2)
	nck = nq - nact # number of constraints for each k = # of unactuated DOFs ( = nunact)
	Bperp = (I - B*B')[nact+1:end,:] # s.t. Bperp*B = 0
	# Number of various constraints - these are used below to set up the jacobian of g.
	ncunact = N * nck
	ncpolytope = size(Cp,1)
	nctransmission = nonlinTransmission ? 1 : 0
	nctotal = ncunact + ncpolytope + nctransmission # this is the TOTAL number of constraints
	dp2 = copy(dp) # No idea why this was getting modified. Storing a copy seems to work.

	eval_g_pieces(k, Δyk, Δykp1, p) = Bperp * Hk(k, Δyk, Δykp1) * (getpt(m, p)[1])
	function eval_g_ret(x)
		pp = x[1:np]
		Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
		gvec = vcat([eval_g_pieces(k, Δyk(k), Δyk(k+1), pp) for k=1:N]...)
		# Polytope constraint on the params
		if ncpolytope > 0
			gvec = [gvec; Cp * pp - dp2] # must be <= 0
		end
		if nctransmission > 0
			gvec = [gvec; gtransmission(m, pp, σomax)]
		end
		if uinfnorm

		end
		return gvec
	end

	# ----------- Constraint Jac ----------------------------
	# Unactuated error: exploit sparsity in the nc*nx matrix. Each constraint depends on Δyk(k), Δyk(k+1), p
	Dgnnz = ncunact * (2*ny + np)
	if ncpolytope > 0
		Dgnnz += prod(size(Cp)) # Assume Cp as a whole is full
	end
	if nctransmission > 0
		Dgnnz += 2 # The +2 at the end is for the transmission constraint
	end

	"Dg jacobian of the constraint function g() for IPOPT."
	function eval_jac_g(x, imode, row::Vector{Int32}, col::Vector{Int32}, value)
		offs = 0
			
		if imode != :Structure
			Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
			p = x[1:np]
			pbb, Tarrr = paramLumped(m, x[1:np])
			τ1, τ2 = Tarrr

			for k=1:N
				# Now assemble the pieces
				dyk = ForwardDiff.jacobian(yy -> eval_g_pieces(k, yy, Δyk(k+1), p), Δyk(k))
				for i=1:nck
					for j=1:ny
						offs += 1
						value[offs] = dyk[i,j]
					end
				end
				dykp1 = ForwardDiff.jacobian(yy -> eval_g_pieces(k, Δyk(k), yy, p), Δyk(k+1))
				for i=1:nck
					for j=1:ny
						offs += 1
						value[offs] = dykp1[i,j]
					end
				end
				dp = ForwardDiff.jacobian(yy -> eval_g_pieces(k, Δyk(k), Δyk(k+1), yy), p)
				for i=1:nck
					for j=1:np
						offs += 1
						value[offs] = dp[i,j]
					end
				end
			end
			if ncpolytope > 0
				# Cp * param <= dp
				for i=1:size(Cp,1)
					for j=1:size(Cp,2)
						offs += 1
						value[offs] = Cp[i,j]
					end
				end
			end
			if nctransmission > 0
				# transmission
				offs += 1
				value[offs] = -σomax/τ1^2 + 4*σomax^3*τ2/(3*τ1^5) # d/dτ1
				offs += 1
				value[offs] = -σomax/(3*τ1^4) # d/dτ2
			end
		else
			for k=1:N
				for i=1:nck
					for j=1:ny
						offs += 1
						row[offs] = (k-1)*nck + i
						col[offs] = np + (k-1)*ny + j
					end
				end
				for i=1:nck
					for j=1:ny
						offs += 1
						row[offs] = (k-1)*nck + i
						col[offs] = np + (k)*ny + j
					end
				end
				for i=1:nck
					for j=1:np
						offs += 1
						row[offs] = (k-1)*nck + i
						col[offs] = j
					end
				end
			end
			if ncpolytope > 0
				# Cp * param <= dp
				for i=1:size(Cp,1)
					for j=1:size(Cp,2)
						offs += 1
						row[offs] = ncunact + i # goes after the unact constraints
						col[offs] = j # hits the elements of param (first part of x)
					end
				end
			end
			if nctransmission > 0
				# transmission
				offs += 1
				row[offs] = ncunact + ncpolytope + 1
				col[offs] = τinds[1]
				offs += 1
				row[offs] = ncunact + ncpolytope + 1
				col[offs] = τinds[2]
			end
		end
	end

	glimsL = -εunact*ones(ncunact)
	glimsU = εunact*ones(ncunact)
	if size(Cp,1) > 0
		glimsL = [glimsL; -100000 * ones(ncpolytope)] # Scott said having non-inf bounds helps IPOPT
		glimsU = [glimsU; zeros(ncpolytope)] # must be <= 0
	end
	if nctransmission > 0
		glimsL = [glimsL; 0.0]
		glimsU = [glimsU; σamax]
	end

	return nctotal, glimsL, glimsU, eval_g_ret, eval_jac_g, Dgnnz, Bperp
end

"Mode=1 => opt, mode=2 ID. Fext(p) or hold constant.

- εunact -- max error to tolerate in the unactuated rows when trying to match passive dynamics.
- plimsL, plimsU -- box constraint for the params. These are used as variable limits in IPOPT.
- σamax -- actuator strain limit. This is used to constrain the transmission coeffs s.t. the actuator displacement is limited to σamax. The form of the actuator constraint depends on bTrCon.
- Cp, dp -- polytope constraint for params. Can pass Cp=ones(0,X) to not include.
"
function optAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, mode::Int, τinds::Array{Int}, R::Tuple, εunact, plimsL, plimsU, σamax, scaleTraj=1.0, trajAct=true, Cp::Matrix=ones(0,1), dp::Vector=ones(0); Fext_pdep::Bool=false, test=false, testTrajReconstruction=false, kwargs...)
	if test
		affineTest(m, opt, traj, param) # this does not need to be here TODO: remove
	end

	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
	np = length(param)
	npt = length(getpt(m, param)[1])

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, scaleTraj; Fext_pdep=Fext_pdep)

	# IPOPT ---------------------------
	# Options on the types of constraints to include
	nonlinTransmission = true # add a transmission constraint in g()? TODO: remove
	uinfnorm = true
	# Transmission limits imposed by actuator
	σomax = norm([yo(k)[1] for k=1:N], Inf)
	Tmin = σomax/σamax
	if !nonlinTransmission
		plimsL[τinds[1]] = Tmin
	end

	"Variables for IPOPT:
	x = [param; Δy] where Δy is the necessary traj modification for passive dynamics matching. 
	For fully-actuated systems, Δy remains 0.
	"
	nx = np + (N+1)*ny # p,Δy
	xlimsL = -1000 * ones(nx)
	xlimsU = 1000 * ones(nx)
	xlimsL[1:np] = plimsL
	xlimsU[1:np] = plimsU
	
	# IPOPT setup using helper functions
	nctotal, glimsL, glimsU, eval_g_ret, eval_jac_g, Dgnnz, Bperp = paramOptConstraint(m, mode, np, ny, nq, δt, Hk, yo, umeas, B, N, Cp, dp, σamax, σomax, εunact, τinds, uinfnorm; nonlinTransmission=nonlinTransmission)
	eval_g(x::Vector, g::Vector) = g .= eval_g_ret(x)
	eval_f, eval_grad_f = paramOptObjQuadratic(m, mode, np, npt, ny, nq, δt, Hk, yo, umeas, B, N, R)
	
	# Create IPOPT problem
	prob = Ipopt.createProblem(
		nx, # Number of variables
		xlimsL, # Variable lower bounds
		xlimsU, # Variable upper bounds
		nctotal, # Number of constraints
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
	trajnew = reconstructTrajFromΔy(m, opt, traj, yo, Hk, B, prob.x[np+1:end], pnew, trajAct)
	unactErr = eval_g_ret(prob.x)

	if testTrajReconstruction
		# Test traj reconstruction:
		Hk, yo, umeas, B, N = paramAffine(m, opt, trajnew, pnew; Fext_pdep=Fext_pdep)
		eval_g_ret2(p) = vcat([Bperp * Hk(k, zeros(ny), zeros(ny)) * (getpt(m, p)[1]) for k=1:N]...)
		display(unactErr')
		display(eval_g_ret2(pnew)')
		error("Tested")
	end

	return pnew, eval_f, trajnew, unactErr, eval_g_ret
end

