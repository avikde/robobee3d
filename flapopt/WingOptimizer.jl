
module WingOptimizer

using ForwardDiff
using Plots; gr()

# only for this model; otherwise pass in these
const ny = 4
const nu = 1

"Plot a list of traj"
function plotTrajs(t, args...)
	traj = args[1]
	N = (length(traj) - ny) ÷ (ny + nu)
	σt = plot(t, traj[1:ny:(N+1)*ny], marker=:auto, ylabel="act disp [m]")
	Ψt = plot(t, traj[2:ny:(N+1)*ny], marker=:auto, ylabel="hinge ang [r]")
	plot(σt, Ψt, layout=(2,1))
	gui()
end

end
