
module WingOptimizer

using ForwardDiff
using Plots; gr()

include("DirTranForm.jl")
include("Wing2DOF.jl")

"Plot a list of traj"
function plotTrajs(t::Vector, ny::Int, nu::Int, params::Vector, args...)
	traj = args[1]
	N = DirTranForm.getN(traj, ny, nu)
	Ny = (N+1)*ny
	# stroke "angle" = T*y[1] / R
	cbar, T = params
	σt = plot(t, traj[1:ny:(N+1)*ny] * T / (Wing2DOF.R/2), marker=:auto, ylabel="stroke ang [r]", title="timestep=$(round(traj[end]; sigdigits=4))ms; c=$(round(cbar; sigdigits=4))mm, T=$(round(T; sigdigits=4))")
	Ψt = plot(t, traj[2:ny:(N+1)*ny], marker=:auto, ylabel="hinge ang [r]")
	ut = plot(t, [traj[Ny+1:nu:Ny+N*nu];NaN], marker=:auto, ylabel="stroke force [mN]")
	plot(σt, Ψt, ut, layout=(3,1))
	gui()
end

function update(traj)
	
end

end
