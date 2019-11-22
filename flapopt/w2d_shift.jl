
include("w2d_paramopt.jl")

"""Hinge phase shift and test opt"""
function shiftAndOpt(ret, minal, shift=0; kwargs...)
	traj1 = copy(ret["traj"])
	function shiftTraj(i)
		comps = traj1[i:ny:(N+1)*ny]
		if shift >= 0
			traj1[i:ny:(N+1)*ny] = [comps[1+shift:end];comps[1:shift]]
		else
			traj1[i:ny:(N+1)*ny] = [comps[end+shift+1:end];comps[1:end+shift]]
		end
	end
	shiftTraj(2)
	shiftTraj(4)
	# after this the traj will likely not be feasible - need to run opt
	retr = opt1(traj1, ret["param"], 1, minal; kwargs...)
	return [retr["param"]; retr["traj"]]
end

function testManyShifts(retdict0, shifts, minal; kwargs...)
	rets = hcat([shiftAndOpt(retdict0, minal, shift; kwargs...) for shift=shifts]...)
	np = length(param0)
	Nshift = size(rets,2)
	paramss = [rets[1:np,i] for i=1:Nshift]
	trajs = [rets[np+1:end,i] for i=1:Nshift]

	pl1 = plotTrajs(m, opt, trajt, paramss, trajs)
	pl2 = plot(shifts, [p[6]/p[2] for p in paramss], marker=:auto, xlabel="shift", ylabel="Tratio")
	plot(pl1..., pl2)
end
