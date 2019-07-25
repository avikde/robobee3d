
module WingOptimizer

using ForwardDiff
using Plots; gr()

include("DirTranForm.jl")
include("Wing2DOF.jl")

"Plot a list of traj"
function plotTrajs(t::Vector, ny::Int, nu::Int, args...)
	traj = args[1]
	N = DirTranForm.getN(traj, ny, nu)
	Ny = (N+1)*ny
	σt = plot(t, traj[1:ny:(N+1)*ny], marker=:auto, ylabel="act disp [m]", title="timestep = $(traj[end])")
	Ψt = plot(t, traj[2:ny:(N+1)*ny], marker=:auto, ylabel="hinge ang [r]")
	ut = plot(t, [traj[Ny+1:nu:Ny+N*nu];NaN], marker=:auto, ylabel="stroke force [N]")
	plot(σt, Ψt, ut, layout=(3,1))
	gui()
end

function update(traj)
	
end

end
