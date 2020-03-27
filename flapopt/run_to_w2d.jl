
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("w2d_model.jl")

m = Wing2DOFModel(
	kbo = [30, 0],
	ma = 6,
	ka = 240,
	Amp = deg2rad.([90, 140]))
ny, nu = cu.dims(m)

function getInitialParams(itype=0)
	# robobee scale
	return 75, [3.2^2,  # cbar2[mm^2] (area/R)^2
		2.6666, # τ1 (from 3333 rad/m, [Jafferis (2016)])
		0.55, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
		2.5, # wΨ [mm]
		0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
		54.4, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
		0.0758 # dt
	]
end
uampl, param0 = getInitialParams()
σamax = 0.3 # [mm] constant? for robobee actuators

include("w2d_paramopt.jl") # for opt
N, trajt, traj0, opt, Φ0 = initTraj(m, param0, 1; uampl=uampl)

# "Cost function components" ------------------

skew(a) = [0 -a[3] a[2];
	a[3] 0 -a[1];
	-a[2] a[1] 0]
wrench(paero, Faero) = [Faero; skew(paero) * Faero]
function wrenchy(y, param, flip=false)
	paero, Jaero, Faero = w2daero(m, y, param)
	if flip
		paero[2] = -paero[2]
		Faero[2] = -Faero[2]
	end
	return wrench(paero, Faero)
end

function wrenchAt(inp, param)
	freq, uamplL, dcL, uamplR, dcR = inp
	thcoeff = 0.1
	Nn = 100
	yL = createInitialTraj(m, opt, Nn, freq, [1e3, 1e2], param, 0; uampl=uamplL, thcoeff=thcoeff, rawtraj=true, verbose=false, dcoffs=dcL) # ny,N array
	Ntot = size(yL, 2)
	
	yR = createInitialTraj(m, opt, Nn, freq, [1e3, 1e2], param, 0; uampl=uamplR, thcoeff=thcoeff, rawtraj=true, verbose=false, dcoffs=dcR)
	
	Nend = 500 # how many to look at
	totalWrench = zeros(Nend, 6)
	# loop through the traj, and compute rcop, F and return vector
	for k=1:Nend
		totalWrench[k,:] = wrenchy(yL[:,Ntot-Nend+k], param) + wrenchy(yR[:,Ntot-Nend+k], param, true)
	end
	return totalWrench
end

tw = wrenchAt([0.16, 75, 0, 75, 0], param0)

pls = [plot(tw[:,c]) for c=1:6]
plot(pls...)
gui()


## ----

# "Objective to minimize"
# function cu.robj(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray)::AbstractArray
#     ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
# 	yk = k -> @view traj[liy[:,k]]
# 	uk = k -> @view traj[liu[:,k]]
    
# 	cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
	
# 	function wrenchk(k)
# 		paero, _, Faero = w2daero(m, cu.transmission(m, yk(k), param)[1], param)
# 		return wrench(paero, Faero)
# 	end

# 	# Compute the average wrench over the cycle
# 	wrenchAvg = zeros(3)
# 	for k=1:N
# 		wrenchAvg += wrenchk(k)
# 	end
# 	wrenchAvg /= N
#     # avg lift
#     return [wrenchAvg[2] - 100]
# end

# # traj opt ------------------------------------

# εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1)
# traj1 = prob.x

# pls = plotTrajs(m, opt, [param0, param0], [traj0, traj1])
# plot(pls...)

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
