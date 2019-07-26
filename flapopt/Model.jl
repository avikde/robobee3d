
using LinearAlgebra, ForwardDiff

"""
Implement these things
"""
abstract type Model end

function dims(m::Model)::Tuple{Int, Int}
	return 0, 0 # ny, nu
end

function dydt(m::Model, y::Vector, u::Vector, params::Vector)::Vector
	return similar(y)
end

function limits(m::Model)::Tuple{Vector, Vector, Vector, Vector}
	ny, nu = dims(m)
	return -Inf*ones(nu), Inf*ones(nu), -Inf*ones(ny), Inf*ones(ny)
end

function limitsTimestep(m::Model)::Tuple{Float64, Float64}
	return 0, Inf
end

function eval_f(m::Model, traj::Vector, params::Vector, N::Int)::Float64
	return 0
end

"Plot a list of traj"
function plotTrajs(m::Model, t::Vector, params::Vector, args...)
	
end

#===========================================================================
Functions that can be specialized optionally
===========================================================================#

"Use autograd to find Jacobians; specialization can do something else"
function Df(m::Model, y::Vector, u::Vector, params::Vector)::Tuple{Matrix, Matrix}
	fy(yy::Vector) = dydt(m, yy, u, params)
	fu(uu::Vector) = dydt(m, y, uu, params)
	Jac_fy = yy -> ForwardDiff.jacobian(fy, yy)
	Jac_fu = uu -> ForwardDiff.jacobian(fu, uu)
	return Jac_fy(y), Jac_fu(u)
end

#===========================================================================
Functions valid for all instances without specialization
===========================================================================#

# Direct transcription helpers
# dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

"Go from traj length"
function Nknot(m::Model, traj::Vector; vart::Bool=true)::Int
	ny, nu = dims(m)
	(length(traj) - ny - (vart ? 1 : 0)) ÷ (ny + nu)
end

"Go from #knot points to traj length"
function Ntraj(m::Model, N::Int; vart::Bool=true)::Int
	ny, nu = dims(m)
	return (N+1)*ny + N*nu + (vart ? 1 : 0)
end

"Return y(k),u(k) for traj, for k ∈ [1,...,(N+1)]"
function linind(m::Model, N::Int)
	ny, nu = dims(m)
	liny = LinearIndices((1:ny, 1:(N+1)))
	# FIXME: the type of this becomes Array due to the +
	linu = LinearIndices((1:nu, 1:(N))) .+ liny[end]
	return liny, linu
end

"Get upper/lower bound on the dirtran variable"
function xbounds(m::Model, N::Int; vart::Bool=true)::Tuple{Vector, Vector}
	umin, umax, xmin, xmax = limits(m)
	x_L = [repeat(xmin, N+1); repeat(umin, N)]
	x_U = [repeat(xmax, N+1); repeat(umax, N)]
	if vart
		δtmin, δtmax = limitsTimestep(m)
		x_L = [x_L; δtmin]
		x_U = [x_U; δtmax]
	end
	return x_L, x_U
end

"""
Dynamics information at state y, input u
	
Note that the actual constraints are:
	g = -ynext + (y + δt * fy)
	dg_dynext = -I
	dg_dy = δt * df_dy + I
	dg_du = δt * df_du
	dg_dδt = fy

Preallocate ouputs like this:
	fy = zeros(ny)
	df_dy = zeros(ny, ny)
	df_du = zeros(ny, nu)
"""
function gdyn!(fy::Vector, df_dy::Matrix, df_du::Matrix, m::Model, y::Vector, u::Vector, params::Vector; order::Int=1)
	fy[:] = dydt(m, y, u, params) # must be constrained to 0
	# Get jacobians wrt VF
	df_dy[:], df_du[:] = Df(m, y, u, params)
end
