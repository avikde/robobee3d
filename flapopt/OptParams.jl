
include("OptBase.jl") #< including this helps vscode reference the functions in there
using Parameters, ForwardDiff, LinearAlgebra, Ipopt, DSP

@with_kw struct ParamOptOpts
	τinds::Array{Int}
	R::Tuple
	plimsL::Vector
	plimsU::Vector
	εunact::Float64 = 0.1
	Fext_pdep::Bool = true
	uinfnorm::Bool = false # only in mode 1
	nonlintransmission::Bool = true # false for linear transmission; true for the cubic polynomial transmission function
end

# --------------
"Implement this"
function paramAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::ParamOptOpts, scaleTraj=1.0)
	error("Implement this!")
end

"Implement this"
paramLumped(m::Model, param::AbstractArray) = error("Implement this")

# lumped parameter vector
function getpt(m::Model, p)
	pb, Tarr, dt = paramLumped(m, p)
	τ1, τ2 = Tarr
	# Nonlinear transmission: see https://github.com/avikde/robobee3d/pull/92
	ptWithTransmission = [pb*τ1; pb*τ2/τ1^2; 1/τ1; τ2/τ1^4]
	# Proper nondim wrt. time  https://github.com/avikde/robobee3d/pull/119#issuecomment-577350049
	return [ptWithTransmission/dt^2; ptWithTransmission/dt; ptWithTransmission], Tarr
end

"Helper function for nonlinear transmission change to H"
Hτ(Hh, σo) = hcat(Hh[:,1:end-2], Hh[:,1:end-2]*σo^2, Hh[:,end-1:end])

"Override this. Input y can be in actuator or output coordinates. If o2a is false, maps A to O, else maps O to A"
function transmission(m::Model, y::AbstractArray, _param::Vector; o2a=false)
	y2 = y
	T = 1.0
	τfun = x -> x
	τifun = x -> x
    return y2, T, τfun, τifun
end

"Some optimization numerical error is causing spiky oscillations in Δy
https://github.com/avikde/robobee3d/pull/127#issuecomment-580050283"
function smoothΔy(Δy, ny, N, dt; ord=2, cutoff_freq=0.5)
	# For each coord, fit a spline
	# tvec = collect(0:(N))*dt
	myfilter = digitalfilter(Lowpass(cutoff_freq; fs=1/dt), Butterworth(ord))
	smoothCoord(i) = filtfilt(myfilter, Δy[i:ny:(N+1)*ny])
	trajNny = hcat([smoothCoord(i) for i=1:ny]...)
	return trajNny'[:]
end

"Helper function to reconstruct the traj (also test it)."
function reconstructTrajFromΔy(m::Model, opt::OptOptions, POPTS, traj::AbstractArray, yo, Δy, paramold, pnew; test::Bool=false)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
	np = length(pnew)
	dtold = paramold[end] # dt is the last param
	dtnew = pnew[end]
	Δy = smoothΔy(Δy, ny, N, dtnew)
	Δyk = k -> Δy[(k-1)*ny+1 : k*ny]

	# Also convert the output traj with the Δy, new T, and inputs
	traj2 = copy(traj)
	# Calculate the new traj (which is in act coordinates, so needs scaling by T)
	ptnew, Tnew = getpt(m, pnew)
	for k=1:N+1
		if opt.trajAct
			# Go from output to act coords
			ya, Tk = transmission(m, yo(k) + Δyk(k), pnew; o2a=true)[1:2]
			yknew = ya
		else
			yknew = yo(k) + Δyk(k)
		end
		# velocity scaling for https://github.com/avikde/robobee3d/issues/110 FIXME: after Δy?
		if abs(dtold/dtnew - 1.0) > 0.001
			yknew[nq+1:end] = yknew[nq+1:end]*dtold/dtnew
		end
		traj2[(k-1)*ny+1:k*ny] = yknew
	end

	# Need to get the affine form again to get the new u. scaleTraj not needed since the yo had it built in
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj2, pnew, POPTS)
	dely0 = zeros(ny)
	for k=1:N
		# Calculate the new inputs
		if k <= N
			traj2[(N+1)*ny+(k-1)*nu+1:(N+1)*ny+(k)*nu] = B' * Hk(k, dely0, dely0)[1] * ptnew
		end
	end

	# println("Test that the lower rows dyn constraint worked ", vcat([Bperp * Hk(k, Δyk(k), Δyk(k+1)) * ptnew for k=1:N]...) - eval_g_ret(prob.x))

	return traj2
end

"""Compute the states at the collocation times tc1 between t1, t2, etc. Note that it uses param provided to evaluate dynamics."""
function collocationStates(m::Model, opt::OptOptions, y, ynext, u, unext, param, dt)
	# precompute state derivatives evaluating the dynamics (USES PARAM)
	dy = dydt(m, y, u, param)
	dynext = dydt(m, ynext, unext, param)

	# Using the formulae in http://underactuated.mit.edu/underactuated.html?chapter=trajopt
	yc = 1/2 * (y + ynext) + dt/8 * (dy - dynext)
	dyc = -3/(2*dt) * (y - ynext) - 1/4 * (dy + dynext)

	return yc, dyc
end

"A more lightweight version of reconstructTrajFromΔy that sets Δy=0"
function getTrajU(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::ParamOptOpts)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, POPTS)
	ptnew, Tnew = getpt(m, param)
	Δy0 = zeros(ny)
	return vcat([B' * Hk(k, Δy0, Δy0)[1] * ptnew for k=1:N]...)
end

function affineTest(m, opt, traj, param, POPTS::ParamOptOpts; fixTraj=false)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
	ptTEST, TTEST = getpt(m, param) # NOTE the actual param values are only needed for the test mode
	dt = param[end]

	traj1 = fixTraj ? fixTrajWithDynConst(m, opt, traj, param) : traj

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj1, param, POPTS)
	Hpb = zeros(nq, N)
	Bu = similar(Hpb)
	for k=1:N
		Hpb[:,k] = Hk(k, zeros(ny), zeros(ny))[1] * ptTEST
		Bu[:,k] = B * umeas(k)[1]
	end
	display(Hpb - Bu)
	# plot of u
	pls = [plot(1:N, [Hpb[i,:]   Bu[i,:]], lw=2)  for i=1:nq]
	plot(pls...)
	gui()
	error("Tested")
end

# TODO: improve this https://github.com/avikde/robobee3d/issues/81
function fixTrajWithDynConst(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	dt = param[end]
	
	# Make a new traj where the dynamics constraint is satisfied exactly
	traj1 = copy(traj)
	yk = k -> @view traj1[liy[:,k]]
	uk = k -> @view traj1[liu[:,k]]
	for k=1:N
		traj1[liy[:,k+1]] = yk(k) + dt * dydt(m, yk(k), uk(k), param)
	end
	return traj1
end

"Smooth ramp function for mechanical power (only positive components) https://math.stackexchange.com/questions/3521169/smooth-approximation-of-ramp-function"
smoothRamp(x; ε=0.1) = x/2 * (1 + x / sqrt(x^2 + ε^2))

"Helper function for optAffine. See optAffine for the def of x"
function paramOptObjective(m::Model, POPTS::ParamOptOpts, mode, np, npt, ny, δt, Hk, yo, umeas, B, N)
	uinfnorm = mode == 2 ? false : POPTS.uinfnorm # no infnorm for ID
	nq = ny÷2
	nact = size(B, 2)
	nΔy = (N+1)*ny
	Δy0 = zeros(ny)
	npt1 = npt÷3 # each of the 3 segments for nondim time https://github.com/avikde/robobee3d/pull/119
	lse = true
	
	Ryy, Ryu, Ruu, wΔy, wu∞ = POPTS.R # NOTE Ryu is just weight on mech. power
	
	function eval_f(x)
		pt, Tarr = getpt(m, x[1:np])
		Δy = x[np+1 : np+nΔy]
		# min Δy
		J = wΔy * dot(Δy, Δy)/N
		if uinfnorm
			s = x[np+nΔy+1 : np+nΔy+nact] # slack variable for infnorm
			J += wu∞ * dot(s, s)
		end
		Jlse = zero(eltype(x))
		
		for k=1:N
			# NOT INCLUDING the Δy in the calculation of u. Including these was resulting in a lot of IPOPT iterations and reconstruction failed -- need to investigate why. https://github.com/avikde/robobee3d/pull/80#issuecomment-541350179
			Hh, Hvel = Hk(k, Δy0, Δy0) #Hk(k, Δyk(k), Δyk(k+1))
			uk = B' * Hh * pt
			if mode == 1
				# For mech pow https://github.com/avikde/robobee3d/issues/123 but use a ramp mapping to get rid of negative power
				# The terms should go in the second segment (/dt) and the last two in that segment (mult by T^-1 terms)
				# ASSUMING 1 ACTUATED DOF
				dqact = [zeros(1,npt1)   Hvel   zeros(1,npt1)] * pt
				J += 1/(2*N) * smoothRamp(dqact' * Ryu * uk)

				# ||u||p
				if !lse
					J += 1/(2*N) * uk' * Ruu * uk
				end
				Jlse += sum(exp.(uk))
			elseif mode == 2
				# ||u||2
				uerr = uk - umeas(k)
				J += 1/(2*N) * uerr' * Ruu * uerr
			end
		end
		if lse
			J += log(Jlse) * Ruu
		end

		return J
	end
	eval_grad_f(x, grad_f) = ForwardDiff.gradient!(grad_f, eval_f, x)

	return eval_f, eval_grad_f
end

"Constraints for IPOPT: g = [gunact; gpolycon; gtransmission]
- gunact = unactuated error = #unactuated DOFS * N
- gpolycon = Cp * p <= dp. But the user can pass in (and default is) Cp 0xX => has no effect
- gtransmission: actuator strain limit (see below)
- ginfnorm: if min of uinfnorm is desired, add on a slack variable s, and add constraints that -s <= uk <= s => {uk-s <= 0, -uk-s <= 0}

Transmission constraint:
	full nonlinear: σomax/τ1 - σomax^3/3 * τ2/τ1^4
	start from a (τ1,0) feasible, then use (see https://github.com/avikde/robobee3d/pull/96#issuecomment-553092480) (τ1-Δτ1, τ2 >= 3/σamax^2*Δτ1) <=> τ2 >= 3/σamax^2*(τ1min - τ1)
	τ2>=0 is in xlims
	The ineq above is: τ1min - τ1 <= τ2*σamax^2/3 => -τ1 - τ2*σamax^2/3 + τ1min <= 0 which is a linear constraint
"
function paramOptConstraint(m::Model, POPTS::ParamOptOpts, mode, np, ny, δt, Hk, yo, umeas, B, N, Cp, dp, σamax, σomax)
	uinfnorm = mode == 2 ? false : POPTS.uinfnorm # no infnorm for ID
	# Unactuated constraint: Bperp' * H(y + Δy) * pt is small enough (unactuated DOFs) 
	nact = size(B, 2)
	nq = ny÷2
	nck = nq - nact # number of constraints for each k = # of unactuated DOFs ( = nunact)
	Bperp = (I - B*B')[nact+1:end,:] # s.t. Bperp*B = 0
	# Number of various constraints - these are used below to set up the jacobian of g.
	ncunact = N * nck
	ncpolytope = size(Cp,1)
	nctransmission = POPTS.nonlintransmission ? 1 : 0
	ncuinfnorm = uinfnorm ? 2 * N * nact : 0
	nctotal = ncunact + ncpolytope + nctransmission + ncuinfnorm # this is the TOTAL number of constraints
	dp2 = copy(dp) # No idea why this was getting modified. Storing a copy seems to work.
	nΔy = (N+1)*ny
	τ1min = σomax/σamax
	# Use delta y in upred (to get uinfnorm to be exact)
	ukpredUseΔy = false
	Δy0 = zeros(ny)

	eval_g_pieces(k, Δyk, Δykp1, p) = Bperp * Hk(k, Δyk, Δykp1)[1] * (getpt(m, p)[1])
	if uinfnorm
		ukpredΔy(k, Δyk, Δykp1, p) = B' * Hk(k, Δyk, Δykp1)[1] * (getpt(m, p)[1])
		ukpred(k, p) = B' * Hk(k, Δy0, Δy0)[1] * (getpt(m, p)[1])
	end

	function eval_g_ret(x)
		pp = x[1:np]
		Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
		gvec = vcat([eval_g_pieces(k, Δyk(k), Δyk(k+1), pp) for k=1:N]...)
		ukp2(k) = ukpredUseΔy ? ukpredΔy(k, Δyk(k), Δyk(k+1), pp) : ukpred(k, pp)
		# Polytope constraint on the params
		if ncpolytope > 0
			gvec = [gvec; Cp * pp - dp2] # must be <= 0
		end
		if nctransmission > 0
			τ1, τ2 = paramLumped(m, pp)[2] # Get both transmission coeffs
			gvec = [gvec; -τ1 - τ2*σamax^2/3 + τ1min]
		end
		if uinfnorm
			s = x[np+nΔy+1 : np+nΔy+nact] # slack variable for infnorm
			gvec = [gvec; 
			vcat([ukp2(k)-s for k=1:N]...); # uk-s<=0
			vcat([-ukp2(k)-s for k=1:N]...)] # -uk-s<=0
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
	if uinfnorm
		# for each k, get 
		# - d/dp(B' * Hk * pt) which should be nact*np, 
		# - identity for the "s", 2x times
		# - d/d(delyk)(B' * Hk * pt) which should be nact*ny * 2 (Δyk, Δykp1)
		Dgnnz += 2 * N * nact * (np + (ukpredUseΔy ? 2 * ny : 0) + 1)
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
				value[offs] = -1 # d/dτ1
				offs += 1
				value[offs] = -σamax^2/3 # d/dτ2
			end
			if uinfnorm
				for k=1:N
					if ukpredUseΔy
						dp = ForwardDiff.jacobian(pp -> ukpredΔy(k, Δyk(k), Δyk(k+1), pp), p)
						dyk = ForwardDiff.jacobian(yy -> ukpredΔy(k, yy, Δyk(k+1), p), Δyk(k))
						dykp1 = ForwardDiff.jacobian(yy -> ukpredΔy(k, Δyk(k), yy, p), Δyk(k+1))
					else
						dp = ForwardDiff.jacobian(pp -> ukpred(k, pp), p)
					end
					for i=1:nact
						for j=1:np
							offs += 1
							value[offs] = dp[i,j] # uk in uk-s<=0
							offs += 1
							value[offs] = -dp[i,j] # uk in -uk-s<=0
						end
						offs += 1
						value[offs] = -1 # -s in uk-s<=0
						offs += 1
						value[offs] = -1 # -s in -uk-s<=0
						if ukpredUseΔy
							for j=1:ny
								offs += 1
								value[offs] = dyk[i,j] # duk/dyk in uk-s<=0
								offs += 1
								value[offs] = -dyk[i,j] # duk/dyk in -uk-s<=0
								offs += 1
								value[offs] = dykp1[i,j] # duk/dykp1 in uk-s<=0
								offs += 1
								value[offs] = -dykp1[i,j] # duk/dykp1 in -uk-s<=0
							end
						end
					end
				end
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
				col[offs] = POPTS.τinds[1]
				offs += 1
				row[offs] = ncunact + ncpolytope + 1
				col[offs] = POPTS.τinds[2]
			end
			if uinfnorm
				rowOffs = ncunact + ncpolytope + nctransmission
				for k=1:N
					for i=1:nact
						for j=1:np
							offs += 1
							row[offs] = rowOffs + (k-1)*nact + i
							col[offs] = j # uk in uk-s<=0
							offs += 1
							row[offs] = rowOffs + ncuinfnorm÷2 + (k-1)*nact + i
							col[offs] = j # uk in -uk-s<=0
						end
						offs += 1
						row[offs] = rowOffs + (k-1)*nact + i
						col[offs] = np+nΔy+i # s in uk-s<=0
						offs += 1
						row[offs] = rowOffs + ncuinfnorm÷2 + (k-1)*nact + i
						col[offs] = np+nΔy+i # s in -uk-s<=0
						if ukpredUseΔy
							for j=1:ny
								offs += 1
								row[offs] = rowOffs + (k-1)*nact + i
								col[offs] = np + (k-1)*ny + j # duk/dyk in uk-s<=0
								offs += 1
								row[offs] = rowOffs + ncuinfnorm÷2 + (k-1)*nact + i
								col[offs] = np + (k-1)*ny + j # duk/dyk in -uk-s<=0
								offs += 1
								row[offs] = rowOffs + (k-1)*nact + i
								col[offs] = np + (k)*ny + j # duk/dykp1 in uk-s<=0
								offs += 1
								row[offs] = rowOffs + ncuinfnorm÷2 + (k-1)*nact + i
								col[offs] = np + (k)*ny + j # duk/dykp1 in -uk-s<=0
							end
						end
					end
				end
			end
		end
	end

	glimsL = -POPTS.εunact*ones(ncunact)
	glimsU = POPTS.εunact*ones(ncunact)
	if size(Cp,1) > 0
		glimsL = [glimsL; -100000 * ones(ncpolytope)] # Scott said having non-inf bounds helps IPOPT
		glimsU = [glimsU; zeros(ncpolytope)] # must be <= 0
	end
	if nctransmission > 0
		glimsL = [glimsL; -100000]
		glimsU = [glimsU; 0.0]
	end
	if uinfnorm
		glimsL = [glimsL; -1e10 * ones(ncuinfnorm)]
		glimsU = [glimsU; zeros(ncuinfnorm)]
	end

	return nctotal, glimsL, glimsU, eval_g_ret, eval_jac_g, Dgnnz, Bperp
end

"Mode=1 => opt, mode=2 ID. Fext(p) or hold constant.

- εunact -- max error to tolerate in the unactuated rows when trying to match passive dynamics.
- plimsL, plimsU -- box constraint for the params. These are used as variable limits in IPOPT.
- σamax -- actuator strain limit. This is used to constrain the transmission coeffs s.t. the actuator displacement is limited to σamax. The form of the actuator constraint depends on bTrCon.
- Cp, dp -- polytope constraint for params. Can pass Cp=ones(0,X) to not include.
"
function optAffine(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::ParamOptOpts, mode::Int, σamax; test=false, testTrajReconstruction=false, Cp::Matrix=ones(0,1), dp::Vector=ones(0), scaleTraj=1.0, dtFix=false, kwargs...)
	if test
		affineTest(m, opt, traj, param, POPTS)
	end
	uinfnorm = mode == 2 ? false : POPTS.uinfnorm # no infnorm for ID
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	np = length(param)
	npt = length(getpt(m, param)[1])
	nΔy = (N+1)*ny

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj, param, POPTS, scaleTraj)
	# IPOPT ---------------------------
	# Options on the types of constraints to include
	σomax = norm([yo(k)[1] for k=1:N], Inf) # Transmission limits imposed by actuator

	"Variables for IPOPT:
	x = [param; Δy] where Δy is the necessary traj modification for passive dynamics matching. 
	For fully-actuated systems, Δy remains 0.
	"
	nx = np + nΔy + (uinfnorm ? nu : 0) # p,Δy[,s]
	xlimsL = -1000 * ones(nx)
	xlimsU = 1000 * ones(nx)
	xlimsL[1:np] = POPTS.plimsL
	xlimsU[1:np] = POPTS.plimsU
	if uinfnorm
		fill!(xlimsL[np+nΔy+1:np+nΔy+nu], 0)
	end
	if dtFix # constrain dt to be the same as it is now
		xlimsL[np] = xlimsU[np] = param[np]
	end
	
	# IPOPT setup using helper functions
	nctotal, glimsL, glimsU, eval_g_ret, eval_jac_g, Dgnnz, Bperp = paramOptConstraint(m, POPTS, mode, np, ny, δt, Hk, yo, umeas, B, N, Cp, dp, σamax, σomax)
	eval_g(x::Vector, g::Vector) = g .= eval_g_ret(x)
	eval_f, eval_grad_f = paramOptObjective(m, POPTS, mode, np, npt, ny, δt, Hk, yo, umeas, B, N)
	
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
	Ipopt.addOption(prob, "sb", "yes") # suppress banner
	for (k,v) in pairs(kwargs)
		# println(k, " => ", v)
		Ipopt.addOption(prob, string(k), v)
	end

	# Solve
	prob.x = [copy(param); zeros(nx - np)]
	status = Ipopt.solveProblem(prob)
	pnew = prob.x[1:np]
	trajnew = reconstructTrajFromΔy(m, opt, POPTS, traj, yo, prob.x[np+1:np+nΔy], param, pnew)

	if testTrajReconstruction
		# Test traj reconstruction:
		Hk, yo, umeas, B, N = paramAffine(m, opt, trajnew, pnew, POPTS)
		Hh = Hk(k, zeros(ny), zeros(ny))[1]
		eval_g_ret2(p) = vcat([Bperp * Hh * (getpt(m, p)[1]) for k=1:N]...)
		display(eval_g_ret2(pnew)')
		error("Tested")
	end
	s = uinfnorm ? prob.x[np+nΔy+1:np+nΔy+nu] : NaN

	return Dict("x"=>prob.x, "traj"=>trajnew, "param"=>prob.x[1:np], "eval_f"=>eval_f, "eval_g"=>eval_g_ret, "s"=>s, "status"=>status)
end

