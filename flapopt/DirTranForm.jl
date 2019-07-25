
module DirTranForm

# dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

# Functions to help with dirtran packing
getN(traj::Vector, ny::Int, nu::Int)::Int = (length(traj) - ny - 1) ÷ (ny + nu)

getSize(N::Int, ny::Int, nu::Int)::Int = (N+1)*ny + N*nu + 1

"Return y(k),u(k) for traj, for k ∈ [1,...,(N+1)]"
function linind(traj::Vector, ny::Int, nu::Int)
	N = getN(traj, ny, nu)
	liny = LinearIndices((1:ny, 1:(N+1)))
	linu = LinearIndices((1:nu, 1:(N))) .+ liny[end]
	return liny, linu
end

"Get upper/lower bound on the dirtran variable"
function getXLU(N::Int, umin, umax, xmin, xmax, δtmin::Float64, δtmax::Float64)
	return [repeat(xmin, N+1); repeat(umin, N); δtmin], [repeat(xmax, N+1); repeat(umax, N); δtmax]
end

end
