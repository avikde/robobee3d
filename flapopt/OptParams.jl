
include("Model.jl") #< including this helps vscode reference the functions in there
using Parameters, ForwardDiff, LinearAlgebra, Ipopt, DSP, SparseArrays

@with_kw mutable struct ParamOptOpts
	τinds::Array{Int}
	"(Ryy, Ryu (mech pow), Ruu, wΔy, wlse)"
	R::Tuple
	plimsL::Vector
	plimsU::Vector
	εunact::Float64 = 0.1
	Fext_pdep::Bool = true
	"Allowed violation for polytope constraint"
	εpoly::Float64 = 1e-4
	"Let the u calculated in the objective depend on Δy -- more accurate but demanding"
	objDepΔy::Bool = true
	"Set positive to enable this constraint (keep diff of successive states in Δy within this bound)."
	ΔySpikyBound::Float64 = 0.0
	"Set to an array of desired params - will use a quadratic weight"
	pdes::Array{Float64} = []
	"Set to an array of (diagonal) weights"
	pdesQ::Array{Float64} = []
	"Use central 1st order difference to avoid phase shift https://github.com/avikde/robobee3d/pull/139"
	centralDiff::Bool = false
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
	pb, τ1, τ2, dt = paramLumped(m, p)
	# Nonlinear transmission: see https://github.com/avikde/robobee3d/pull/92
	ptWithTransmission = [pb*τ1; pb*τ2/τ1^2; 1/τ1; τ2/τ1^4]
	# Proper nondim wrt. time  https://github.com/avikde/robobee3d/pull/119#issuecomment-577350049
	return [ptWithTransmission/dt^2; ptWithTransmission/dt; ptWithTransmission]
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
	# Calculate the new traj
	ptnew = getpt(m, pnew)
	for k=1:N+1
		yknew = yo(k) + Δyk(k)
		# velocity scaling for https://github.com/avikde/robobee3d/issues/110 after Δy since that is how Δy was calculated
		if abs(dtold/dtnew - 1.0) > 0.001
			yknew[nq+1:end] = yknew[nq+1:end]*dtold/dtnew
		end
		traj2[(k-1)*ny+1:k*ny] = yknew
	end

	# Need to get the affine form again to get the new u. scaleTraj not needed since the yo had it built in
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj2, pnew, POPTS)
	Δy0 = zeros(ny)
	for k=1:N
		# Calculate the new inputs
		if k <= N
			traj2[(N+1)*ny+(k-1)*nu+1:(N+1)*ny+(k)*nu] = B' * Hk(k, Δy0, Δy0, Δy0)[1] * ptnew
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
	ptnew = getpt(m, param)
	Δy0 = zeros(ny)
	return vcat([B' * Hk(k, Δy0, Δy0, Δy0)[1] * ptnew for k=1:N]...)
end

function affineTest(m, opt, traj, param, POPTS::ParamOptOpts; fixTraj=false)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	nq = ny÷2
	ptTEST = getpt(m, param) # NOTE the actual param values are only needed for the test mode
	dt = param[end]
	Δy0 = zeros(ny)

	traj1 = fixTraj ? fixTrajWithDynConst(m, opt, traj, param) : traj

	# Quadratic form matrix
	Hk, yo, umeas, B, N = paramAffine(m, opt, traj1, param, POPTS)
	Hpb = zeros(nq, N)
	Bu = similar(Hpb)
	for k=1:N
		Hpb[:,k] = Hk(k, Δy0, Δy0, Δy0)[1] * ptTEST
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
ramp(x; ε=0.1, smooth=true) = smooth ? (x/2 * (1 + x / sqrt(x^2 + ε^2))) : (x > 0 ? x : 0)
dramp(x; ε=0.1, smooth=true) = smooth ? (1/2*(1 + (x*(x^2 + 2*ε^2))/(x^2 + ε^2)^(3/2))) #= used Mathematica =# : sign(x) 
"https://en.wikipedia.org/wiki/LogSumExp"
LSE(x) = log(sum(exp.(x)))
function dLSE(x)
	y = exp.(x)
	return y/sum(y)
end

"Assemble the big Hu,Hdq matrices s.t. Hu * pt = uact, Hdq * pt = dqact - ASSUMING dely = 0.
NOT INCLUDING the Δy in the calculation of u. Including these was resulting in a lot of IPOPT iterations and reconstruction failed -- need to investigate why. https://github.com/avikde/robobee3d/pull/80#issuecomment-541350179"
function bigH(N, ny, nact, npt, Hk, B, Δy)
	npt1 = npt÷3 # each of the 3 segments for nondim time https://github.com/avikde/robobee3d/pull/119
	Δyk = k -> Δy[(k-1)*ny+1 : k*ny]
	nq = ny÷2
	Hu = zeros(nact*N, npt)
	Hdq = similar(Hu)
	# Bperp = (I - B*B')[nact+1:end,:] # s.t. Bperp*B = 0
	# nunact = nq - nact
	# Hunact = zeros(nunact*N, npt)
	for k=1:N
		Hh, Hvel = Hk(k, Δyk(max(k-1,1)), Δyk(k), Δyk(k+1))
		Hu[nact*(k-1)+1 : nact*k, :] = B' * Hh
		# Hunact[nunact*(k-1)+1 : nunact*k, :] = Bperp * Hh
		# The terms should go in the second segment (/dt) and the last two in that segment (mult by T^-1 terms)
		# ASSUMING 1 ACTUATED DOF
		Hdq[nact*(k-1)+1 : nact*k, :] = [zeros(1,npt1)  Hvel  zeros(1,npt1)]
	end
	return Hu, Hdq#, Hunact
end

"Helper function for paramOpt. See paramOpt for the def of x"
function paramOptObjective(m::Model, POPTS::ParamOptOpts, mode, np, npt, ny, δt, Hk, yo, umeas, B, N)
	dφ_dptAutodiff = true
	nq = ny÷2
	nact = size(B, 2)
	nΔy = (N+1)*ny
	objDepΔy = nact == nq ? false : POPTS.objDepΔy # If fully actuated no Δy
	
	Ryy, Ryu, Ruu, wΔy, wlse = POPTS.R # NOTE Ryu is just weight on mech. power
	lse = wlse > 1e-6
	NJcomps = 5
	usePdes = length(POPTS.pdes) > 0

	if !objDepΔy # Ignore Δy in computation of H for objective (exactly fine for fully act)
		Hu, Hdq = bigH(N, ny, nact, npt, Hk, B, zeros(nΔy))
	end

	"Running cost function of actuator force, vel for a single k"
	function φmech(uact, dqact)
		Jcomps = zeros(eltype(uact), NJcomps)
		if mode == 1
			# For mech pow https://github.com/avikde/robobee3d/issues/123 but use a ramp mapping to get rid of negative power
			Jcomps[3] = 1/2 * ramp(dqact' * Ryu * uact)

			# ||u||p
			Jcomps[4] = 1/2 * uact' * Ruu * uact
		elseif mode == 2
			# ||u||2
			uerr = uact - umeas(k)
			Jcomps[4] = 1/2 * uerr' * Ruu * uerr
		end
		return Jcomps
	end
	"Analytical gradient of the above (of the summed cost). No LSE yet"
	function dφmech(uact, dqact)
		if mode == 1
			dsmooth = 1/2 * dramp(dqact' * Ryu * uact)
			return vcat(Ruu * uact + dsmooth * Ryu * dqact, dsmooth * Ryu * uact)
		elseif mode == 2
			return vcat(Ruu * (uact - umeas(k)), zero(uact))
		end
	end

	unpackX(x) = getpt(m, x[1:np]), x[np+1 : np+nΔy]

	_k(k) = nact*(k-1)+1 : nact*k

	function φ(p, pt, Δy; debug=false)
		# min Δy
		T = eltype(pt)
		Jcomps = zeros(T, NJcomps) # [JΔy, Jlse, Jpow, Ju2]

		if objDepΔy
			Hu, Hdq = bigH(N, ny, nact, npt, Hk, B, Δy) # use the new Δy
		end
		uvec = Hu * pt
		dqvec = Hdq * pt

		for k=1:N
			Jcomps .+= φmech(uvec[_k(k)], dqvec[_k(k)])
		end

		# Total
		Jcomps[1] = wΔy/N * dot(Δy, Δy)
		if lse
			Jcomps[2] = wlse * LSE(uvec)
		end
		Jcomps[3] /= N
		Jcomps[4] /= N
		if usePdes
			perr = p - POPTS.pdes
			Jcomps[5] = 1/2 * perr' * Diagonal(POPTS.pdesQ) * perr
		end
		if debug
			return pt, Hk, B, Jcomps, [uvec  dqvec]
		end

		return sum(Jcomps)
	end
	eval_f(x; kwargs...) = φ(x[1:np], unpackX(x)...; kwargs...)

	#Storage for ForwardDiff
	dφ_dpt = zeros(npt)
	dpt_dp = zeros(npt, np)

	function eval_grad_f(x, grad_f)
		pt, Δy = unpackX(x)
		ForwardDiff.jacobian!(dpt_dp, x -> getpt(m, x), x[1:np]) # 1,npt X npt,np

		if dφ_dptAutodiff
			# Autodiff for gradient of 
			ForwardDiff.gradient!(dφ_dpt, ptdiff -> φ(x[1:np], ptdiff, Δy), pt)
			grad_f[1:np] = dφ_dpt' * dpt_dp
		else
			# Analytical gradients https://github.com/avikde/robobee3d/pull/137
			if objDepΔy
				Hu, Hdq = bigH(N, ny, nact, npt, Hk, B, Δy)
			end
			uvec = Hu * pt
			dqvec = Hdq * pt
			# Need this first to "clear" previous grad_f (or could fill it with 0)
			grad_f[1:np] = wlse * (dLSE(uvec)' * Hu * dpt_dp) # For LSE
			for k=1:N
				dφmech1 = dφmech(uvec[_k(k)], dqvec[_k(k)])
				grad_f[1:np] += 1/N * (dφmech1' * [Hu[_k(k),:]; Hdq[_k(k),:]] * dpt_dp)[:]
			end
		end

		grad_f[np+1:np+nΔy] = 2*wΔy/N*Δy # nΔy Analytical - see cost above
		if usePdes
			grad_f[1:np] += POPTS.pdesQ .* (x[1:np] - POPTS.pdes) # gradient of quadratic
		end
	end
	
	return eval_f, eval_grad_f
end

"Construct a N*ny x (N+1)*ny matrix D s.t. D*y is the traj diff where each element is yi[k+1]-yi[k]"
function trajDiffMat(N, ny)
	Diffm = spdiagm(0 => ones((N+1)*ny), 4 => -ones(N*ny))
	return Diffm[1:(N*ny), :]
end

"Return in the form C, d s.t. C p <= d.
See Mathematica, linearization gives -τ1 - τ2*σamax^2/3 + τ1min <= 0.
Transmission constraint:
	full nonlinear: σomax/τ1 - σomax^3/3 * τ2/τ1^4
	start from a (τ1,0) feasible, then use (see https://github.com/avikde/robobee3d/pull/96#issuecomment-553092480) (τ1-Δτ1, τ2 >= 3/σamax^2*Δτ1) <=> τ2 >= 3/σamax^2*(τ1min - τ1)
	τ2>=0 is in xlims
	The ineq above is: τ1min - τ1 <= τ2*σamax^2/3 => -τ1 - τ2*σamax^2/3 + τ1min <= 0 which is a linear constraint"
function transmissionLinearConstraint(POPTS, np, σamax, σomax)
	τ1min = σomax/σamax # this would be the case with the linear transmission
	C = zeros(1,np)
	for i=1:np
		if i == POPTS.τinds[1]
			C[i] = -1.0
		elseif i == POPTS.τinds[2]
			C[i] = -σamax^2/3
		end
	end
	return C, [-τ1min]
end

"Constraints for IPOPT: g = [gunact; gpolycon; gtransmission]
- gunact = unactuated error = #unactuated DOFS * N
- gpolycon = Cp * p <= dp. But the user can pass in (and default is) Cp 0xX => has no effect
- gtransmission: actuator strain limit (see below)
"
function paramOptConstraint(m::Model, POPTS::ParamOptOpts, mode, np, ny, δt, Hk, yo, umeas, B, N, Cp, dp)
	# Unactuated constraint: Bperp' * H(y + Δy) * pt is small enough (unactuated DOFs) 
	nact = size(B, 2)
	nq = ny÷2
	nunact = nq - nact # number of constraints for each k = # of unactuated DOFs ( = nunact)
	Bperp = (I - B*B')[nact+1:end,:] # s.t. Bperp*B = 0
	# Number of various constraints - these are used below to set up the jacobian of g.
	ncunact = N * nunact
	
	# Convert the polytope constraint to sparse, and get the COO format IPOPT needs
	CpS = sparse(Cp)
	CpI, CpJ, CpV = findnz(CpS)
	ncpolytope = CpS.m
	ncpolytopennz = length(CpV)

	# Δy diff constraint (to try and avoid spikiness) https://github.com/avikde/robobee3d/pull/138. It works but adds a 
	ΔySpikyConst = POPTS.ΔySpikyBound > 1e-4
	if ΔySpikyConst
		spikyD = trajDiffMat(N, ny)
		DI, DJ, DV = findnz(spikyD)
		ncspiky = spikyD.m
		ncspikynnz = length(DI)
	else
		ncspiky = ncspikynnz = 0
	end

	nctotal = ncunact + ncpolytope + ncspiky # this is the TOTAL number of constraints
	dp2 = copy(dp) # No idea why this was getting modified. Storing a copy seems to work.
	nΔy = (N+1)*ny

	function eval_g(x, g)
		pp = x[1:np]
		pt = getpt(m,pp)
		Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
		
		g .= [vcat([Bperp * Hk(k, Δyk(max(k-1,1)), Δyk(k), Δyk(k+1))[1] * pt for k=1:N]...); # unact
			CpS * pp; # Polytope constraint on the params must be <= dp
			ΔySpikyConst ? spikyD * x[np+1 : np+(N+1)*ny] : []] # D*Δy within bounds
	end

	# ----------- Constraint Jac ----------------------------
	# Unactuated error: exploit sparsity in the nc*nx matrix. Each constraint depends on Δyk(k-1), Δyk(k), Δyk(k+1), p
	Dgnnz = ncunact * (3*ny + np) + ncpolytopennz + ncspikynnz
	# Storage for Jacobians
	dykm1 = zeros(nunact, ny)
	dyk = similar(dykm1)
	dykp1 = similar(dykm1)
	dp = zeros(nunact, np)

	"Dg jacobian of the constraint function g() for IPOPT."
	function eval_jac_g(x, imode, row::Vector{Int32}, col::Vector{Int32}, value)
		offs = 0
		
		if imode == :Values
			Δyk = k -> x[np+(k-1)*ny+1 : np+k*ny]
			p = x[1:np]
			pt = getpt(m, p)
			# can replace by inplace
			dpt_dp = ForwardDiff.jacobian(x -> getpt(m, x), p) # 1,npt X npt,np
		end

		for k=1:N
			# Use AD for these. These individual terms should not be expected to be sparse since it is a np -> nunact map.
			kprev = max(k-1,1)
			if imode == :Values
				# Tested https://github.com/avikde/robobee3d/pull/139 that explicitly leaving this out if not centralDiff (and reducing nnz) does not perform better and is slower?
				ForwardDiff.jacobian!(dykm1, yy -> Bperp * Hk(k, yy, Δyk(k), Δyk(k+1))[1] * pt, Δyk(kprev))
				ForwardDiff.jacobian!(dyk, yy -> Bperp * Hk(k, Δyk(kprev), yy, Δyk(k+1))[1] * pt, Δyk(k))
				ForwardDiff.jacobian!(dykp1, yy -> Bperp * Hk(k, Δyk(kprev), Δyk(k), yy)[1] * pt, Δyk(k+1))
				dp .= Bperp * Hk(k, Δyk(kprev), Δyk(k), Δyk(k+1))[1] * dpt_dp
			end
			for i=1:nunact
				for j=1:ny
					offs += 1
					if imode == :Values
						value[offs] = dykm1[i,j]
					else
						row[offs] = (k-1)*nunact + i
						col[offs] = np + (kprev-1)*ny + j
					end
				end
			end
			for i=1:nunact
				for j=1:ny
					offs += 1
					if imode == :Values
						value[offs] = dyk[i,j]
					else
						row[offs] = (k-1)*nunact + i
						col[offs] = np + (k-1)*ny + j
					end
				end
			end
			for i=1:nunact
				for j=1:ny
					offs += 1
					if imode == :Values
						value[offs] = dykp1[i,j]
					else
						row[offs] = (k-1)*nunact + i
						col[offs] = np + (k)*ny + j
					end
				end
			end
			for i=1:nunact
				for j=1:np
					offs += 1
					if imode == :Values
						value[offs] = dp[i,j]
					else
						row[offs] = (k-1)*nunact + i
						col[offs] = j
					end
				end
			end
		end

		# Cp * param <= dp
		for i=1:ncpolytopennz
			offs += 1
			if imode == :Values
				value[offs] = CpV[i]
			else
				row[offs] = ncunact + CpI[i] # goes after the unact constraints
				col[offs] = CpJ[i] # hits the elements of param (first part of x)
			end
		end

		# D * Δy within bounds
		for i=1:ncspikynnz
			offs += 1
			if imode == :Values
				value[offs] = DV[i]
			else
				row[offs] = ncunact + ncpolytope + DI[i] # goes after the previous constraints
				col[offs] = np + DJ[i] # hits Δy (after param)
			end
		end
	end

	# Append to lims
	glimsL = -POPTS.εunact*ones(ncunact)
	glimsU = POPTS.εunact*ones(ncunact)

	glimsL = [glimsL; -1e3 * ones(ncpolytope)] # Scott said having non-inf bounds helps IPOPT
	glimsU = [glimsU; dp2 + POPTS.εpoly * ones(ncpolytope)] # must be <= 0

	glimsL = [glimsL; -POPTS.ΔySpikyBound * ones(ncspiky)] # empirical from https://github.com/avikde/robobee3d/pull/138
	glimsU = [glimsU; POPTS.ΔySpikyBound * ones(ncspiky)]
	
	return nctotal, glimsL, glimsU, eval_g, eval_jac_g, Dgnnz, Bperp
end

"Mode=1 => opt, mode=2 ID. Fext(p) or hold constant.

- εunact -- max error to tolerate in the unactuated rows when trying to match passive dynamics.
- plimsL, plimsU -- box constraint for the params. These are used as variable limits in IPOPT.
- σamax -- actuator strain limit. This is used to constrain the transmission coeffs s.t. the actuator displacement is limited to σamax. The form of the actuator constraint depends on bTrCon.
- Cp, dp -- polytope constraint for params. Can pass Cp=ones(0,X) to not include.
"
function paramOpt(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::ParamOptOpts, mode::Int, σamax; test=false, Cp::Matrix=ones(0,1), dp::Vector=ones(0), scaleTraj=1.0, dtFix=false, kwargs...)
	if test
		affineTest(m, opt, traj, param, POPTS)
	end
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	np = length(param)
	npt = length(getpt(m, param))
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
	nx = np + nΔy # p,Δy
	xlimsL = -1e3 * ones(nx)
	xlimsU = 1e3 * ones(nx)
	xlimsL[1:np] = POPTS.plimsL
	xlimsU[1:np] = POPTS.plimsU
	if dtFix # constrain dt to be the same as it is now
		xlimsL[np] = xlimsU[np] = param[np]
	end
	
	# IPOPT setup using helper functions
	Ct, dt = transmissionLinearConstraint(POPTS, np, σamax, σomax) # Append linear transmission constraint
	nctotal, glimsL, glimsU, eval_g, eval_jac_g, Dgnnz, Bperp = paramOptConstraint(m, POPTS, mode, np, ny, δt, Hk, yo, umeas, B, N, [Cp; Ct], [dp; dt])
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
	Ipopt.addOption(prob, "sb", "yes") # suppress banner
	Ipopt.addOption(prob, "print_level", 1) # no printing (user can override)

	# Add options using kwargs
	for (k,v) in pairs(kwargs)
		# println(k, " => ", v)
		Ipopt.addOption(prob, string(k), v)
	end

	# Solve
	prob.x = [copy(param); zeros(nx - np)]
	status = Ipopt.solveProblem(prob)
	pnew = prob.x[1:np]
	trajnew = reconstructTrajFromΔy(m, opt, POPTS, traj, yo, prob.x[np+1:np+nΔy], param, pnew)

	return Dict("x"=>prob.x, "traj"=>trajnew, "param"=>prob.x[1:np], "eval_f"=>eval_f, "eval_grad_f"=>eval_grad_f, "eval_g"=>eval_g, "nc"=>nctotal, "status"=>status)
end

