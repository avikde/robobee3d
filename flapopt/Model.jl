

"""
Implement these things
"""
abstract type Model end

@enum OptBoundaryConstraint NONE SYMMETRIC PERIODIC

"An immutable struct of options"
struct OptOptions
	vart::Bool
	fixedδt::Float64 # irrelevant if vart=true
	order::Int
	boundaryConstraint::OptBoundaryConstraint
end

#=========================================================================
Functions that must be specialized by a Model
=========================================================================#

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

function Jobj(m::Model, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3)::Number
	return 0
end

#=========================================================================
Functions that can be specialized optionally
=========================================================================#

"Discrete linearization using autograd (model must provide dydt).
A model can specialize this function to m::MyModel if it is already linear"
function dlin!(Ak::Matrix, Bk::Matrix, m::Model, y::Vector, u::Vector, params::Vector, δt::Float64)
	# Autograd linearization of dydt
	fy(yy::Vector) = dydt(m, yy, u, params)
	fu(uu::Vector) = dydt(m, y, uu, params)
	ForwardDiff.jacobian!(Ak, fy, y)
	ForwardDiff.jacobian!(Bk, fu, u)
	# Continuous -> discrete
	Ak .= I + δt * Ak
	Bk .= δt * Bk
	# TODO: affine term??
end

"Discrete dynamics step"
function ddynamics(m::Model, y::Vector, u::Vector, params::Vector, δt::Float64; useLinearization::Bool=false)
	return y + δt * dydt(m, y, u, params)
	# TODO: useLinearization
end

"Discrete linearization wrt params."
function dlinp!(Pk::Matrix, m::Model, y::Vector, u::Vector, params::Vector, δt::Float64)
	ForwardDiff.jacobian!(Pk, pp -> dydt(m, y, u, pp), params)
	Pk .= δt * Pk
end

#=========================================================================
Functions valid for all instances without specialization
=========================================================================#

# Direct transcription helpers
# dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

"Go from traj length"
function Nknot(m::Model, opt::OptOptions, traj::Vector)::Int
	ny, nu = dims(m)
	(length(traj) - ny - (opt.vart ? 1 : 0)) ÷ (ny + nu)
end

"Go from #knot points to traj length"
function Ntraj(m::Model, opt::OptOptions, N::Int)::Int
	ny, nu = dims(m)
	return (N+1)*ny + N*nu + (opt.vart ? 1 : 0)
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
function xbounds(m::Model, opt::OptOptions, N::Int)::Tuple{Vector, Vector}
	umin, umax, xmin, xmax = limits(m)
	x_L = [repeat(xmin, N+1); repeat(umin, N)]
	x_U = [repeat(xmax, N+1); repeat(umax, N)]
	if opt.vart
		δtmin, δtmax = limitsTimestep(m)
		x_L = [x_L; δtmin]
		x_U = [x_U; δtmax]
	end
	return x_L, x_U
end

"Return the dynamics timestep"
function getδt(opt::OptOptions, traj::Vector)::Number
	return opt.vart ? traj[end] : opt.fixedδt
end

"Helper to get all model info in one line"
function modelInfo(m::Model, opt::OptOptions, traj::Vector)::Tuple
	ny, nu = dims(m)
	N = Nknot(m, opt, traj)
	δt = getδt(opt, traj)
	liy, liu = linind(m, N)
	return ny, nu, N, δt, liy, liu
end
