
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("Wing2DOF.jl")

# create an instance
# From Patrick 300 mN-mm/rad. 1 rad => R/2 σ-displacement. The torque is applied with a lever arm of R/2 => force = torque / (R/2)
# so overall, get 300 / (R/2)^2.
# Noah said lumped stiffness is ~60% transmission and 40% actuator => k = 1.5 = ko + ka/T^2.
# For T=20, get ka = 240.
# To get ma, use the fact that actuator resonance is ~1KHz => equivalent ma = 240/(2*pi)^2 ~= 6mg
m = Wing2DOFModel(
	17.0, # R, [Jafferis (2016)]
	0.55#= 1.5 =#, #k output
	0, #b output
	6, # ma
	0, # ba
	150#= 0 =#, # ka
	true) # bCoriolis
ny, nu = cu.dims(m)

function getInitialParams(itype=0)
	if itype==0
		# robobee scale
		return 70, [3.2,  # cbar[mm] (area/R)
			28.33, # τ1 (from 3333 rad/m, R=17, [Jafferis (2016)])
			0.52, # mwing[mg]
			5, # kΨ [mN-mm/rad]
			3, # bΨ [mN-mm/(rad/ms)]
			0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
			0.135 # dt
		]
	elseif itype==1
		# scaled up
		return 100, [5.411,  # cbar[mm] (area/R)
			18.681, # τ1 (from 3333 rad/m, R=17, [Jafferis (2016)])
			0.866, # mwing[mg]
			22.633, # kΨ [mN-mm/rad]
			12.037, # bΨ [mN-mm/(rad/ms)]
			0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
			0.109 # dt
		]
	end
end
uampl, param0 = getInitialParams(1)
σamax = 0.3 # [mm] constant? for robobee actuators

includet("w2d_paramopt.jl")

# IMPORTANT - load which traj here!!!
KINTYPE = 1
N, trajt, traj0, opt, avgLift0 = initTraj(KINTYPE; uampl=uampl)

# # traj opt ------------------------------------

# # εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# # prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj)

# # plot(plot(prob.g), plot(prob.mult_g), size=(900,400))

# mo = nothing#cu.paramoptQPSetup(m, opt, traj0; scaling=false, verbose=false)

# # IPOPT
# εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1, nlp_scaling_method="none")
# traj1 = prob.x
# # trajs = [traj0, traj1]
# # params = [param0, param0]

# # # naive param opt
# # εsp = [100.0, 0.005, 100.0] # IC, dyn, symm
# # param1 = cu.ipoptsolve(m, opt, traj0, param0, εs, :param)
# # trajs = [traj0, traj0]
# # params = [param0, param1]

# # with my modification to Ha/Coros g-preferred param opt
# δx = cu.paramδx(m, opt, traj0, param0, prob.mult_x_L, prob.mult_x_U)
# param1 = cu.paramopt(mo, m, opt, traj1, param0, δx, εs, Q; step=1e2)
# # param1 = cu.paramoptJ(m, opt, traj1, param0, εs; step=0.01)
# prob = cu.ipoptsolve(m, opt, traj1, param1, εs, :traj; print_level=1, nlp_scaling_method="none")
# traj2 = prob.x

# trajs = [traj0, traj1, traj2]
# params = [param0, param0, param1]

# # # Custom solver ---
# # wkt = cu.OptWorkspace(cu.Ntraj(m, opt, N), (N+2)*ny)
# # @time traj1 = cu.csSolve!(wkt, m, opt, traj0, param0, :traj; Ninner=30, μs=[1e5])
# # wkp = cu.OptWorkspace(length(param0), (N+2)*ny)
# # @time param1 = cu.csSolve!(wkp, m, opt, traj1, param0, :param; Ninner=30, μs=[1e3])
# # @time traj2 = cu.csSolve!(wkt, m, opt, traj1, param1, :traj; Ninner=30, μs=[1e5])
# # trajs = [traj0, traj1, traj2]
# # params = [param0, param0, param1]
# # # trajs, params, wkt = cu.csAlternateSolve(m, opt, traj0, params0, 1; μst=[1e6], Ninnert=30, μsp=[1e-2,1e-2], Ninnerp=2)

# # # pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# # # display(params)

# # paramopt -------------------------------------


# # visualize --------------------------------------------

# println("Objectives: ", [(_ro = cu.robj(m, opt, tt, param0); _ro ⋅ _ro) for tt in trajs])
# println("Params: ", params)

# # animateTrajs(m, opt, params, trajs)
# pl1 = plotTrajs(m, opt, params, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# # plot(pl1...)
# # gui()
