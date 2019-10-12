
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
function gvalues!(gout::AbstractArray, m::Model, opt::OptOptions, traj::AbstractArray, params::AbstractArray, y0::AbstractArray)
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
			gout[li[:,k+1]], fkp1 = gkDirCol(m, params, yk(k), yk(k+1), uk(k), ukp1, δt, fk)
		end
	end

	# Periodicity or symmetry
	if opt.boundaryConstraint == :symmetric
		# FIXME: for now hardcoded symmetry G(y) = -y
		gout[li[:,N+2]] = -yk(1) - yk(N+1)
	end

	return
end

# TODO: improve this https://github.com/avikde/robobee3d/issues/81
function fixTrajWithDynConst(m::Model, opt::OptOptions, traj::AbstractArray, param::AbstractArray)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	
	# Make a new traj where the dynamics constraint is satisfied exactly
	traj1 = copy(traj)
	yk = k -> @view traj1[liy[:,k]]
	uk = k -> @view traj1[liu[:,k]]
	for k=1:N
		traj1[liy[:,k+1]] = yk(k) + δt * dydt(m, yk(k), uk(k), param)
	end
	return traj1
end

#=========================================================================
Visualization
=========================================================================#

function visualizeConstraintViolations(m::Model, opt::OptOptions, params, trajs)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, trajs[1])
	g_L, g_U = gbounds(m, opt, trajs[1])

	function getg(traj, param)
		g0 = similar(g_L)
		gvalues!(g0, m, opt, traj, param, trajs[1][1:ny])
		return g0
	end

	# println(g0, g1)
    Nt = length(trajs)
	pl2 = plot([getg(trajs[i], params[i]) for i=1:Nt], marker=:auto, title="Constraint violations")
	hline!(pl2, [0], color="black", alpha=0.3)
end

