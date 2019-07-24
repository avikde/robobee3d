
module WingOptimizer

using ForwardDiff
using Plots; gr()

# only for this model; otherwise pass in these
const ny = 4
const nu = 1

# Functions to help with dirtran packing
getN = traj -> (length(traj) - ny) ÷ (ny + nu)
"Return y(k),u(k) for traj, for k ∈ [1,...,(N+1)]"
function yuk(traj, k)
	N = getN(traj)
	liny = LinearIndices((1:ny, 1:(N+1)))
	linu = liny[end] .+ LinearIndices((1:nu, 1:(N)))
	return traj[liny[:,k]], traj[linu[:,k]]
end

"Plot a list of traj"
function plotTrajs(t, args...)
	traj = args[1]
	N = getN(traj)
	σt = plot(t, traj[1:ny:(N+1)*ny], marker=:auto, ylabel="act disp [m]")
	Ψt = plot(t, traj[2:ny:(N+1)*ny], marker=:auto, ylabel="hinge ang [r]")
	plot(σt, Ψt, layout=(2,1))
	gui()
end

function update(traj)
	
end

end
