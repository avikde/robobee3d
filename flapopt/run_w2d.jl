
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("Wing2DOF.jl")
includet("LoadWingKinData.jl")

# create an instance
# From Patrick 300 mN-mm/rad. 1 rad => R/2 σ-displacement. The torque is applied with a lever arm of R/2 => force = torque / (R/2)
# so overall, get 300 / (R/2)^2.
# Noah said lumped stiffness is ~60% transmission and 40% actuator => k = 1.5 = ko + ka/T^2.
# For T=20, get ka = 240.
# To get ma, use the fact that actuator resonance is ~1KHz => equivalent ma = 240/(2*pi)^2 ~= 6mg
m = Wing2DOFModel(
	17.0, # R, [Jafferis (2016)]
	0.9#= 1.5 =#, #k output
	0, #b output
	6, # ma
	0, # ba
	150#= 0 =#) # ka
ny, nu = cu.dims(m)
param0 = [3.2,  # cbar[mm] (area/R)
	28.33, # T (from 3333 rad/m, R=17, [Jafferis (2016)])
	0.52, # mwing[mg]
	5, # kΨ [mN-mm/rad]
	3 # bΨ [mN-mm/(rad/ms)]
]

# Stiffness sweep ---
# function respkσ(kσ)
# 	olrfun = f -> openloopResponse(m, opt, f, [param0; kσ])
# 	freqs = collect(0.05:0.01:0.45)
# 	mags = hcat(olrfun.(freqs)...)'
# 	return plot(freqs, mags, ylabel=string(kσ * 100), legend=false)
# end
# kσs = collect(0.0:0.5:5)
# pls = respkσ.(kσs)
# plot(pls...)

# Sim data
opt = cu.OptOptions(false, 0.2, 1, :none, 1e-8, false) # sim
N = opt.boundaryConstraint == :symmetric ? 17 : 33
trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0, 187)

# Load data
# opt = cu.OptOptions(false, 0.1, 1, :none, 1e-8, false) # real
# N, trajt, traj0, lift, drag = loadAlignedData("data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.2445; strokeMult=m.R/(2*param0[2]), ForcePerVolt=0.8)
# pl1 = compareTrajToDAQ(m, opt, trajt, param0, traj0, lift, drag)
# plot(pl1...)

# Make traj satisfy dyn constraint with these params?
traj0 = cu.fixTrajWithDynConst(m, opt, traj0, param0)

# Constraint on cbar placed by minAvgLift. FIXME: this is very specific to W2D, since lift \proptp cbar
avgLift0 = avgLift(m, opt, traj0, param0)
cbarmin = minAvgLift -> param0[1] * minAvgLift / avgLift0

R_WTS = (zeros(4,4), 0, 1.0*I)#diagm(0=>[0.1,100]))

# FIXME: σomax is printed by optAffine
# σamax = 0.3 # [mm] constant? for robobee actuators
Tmin = 19.011058431792932 # σomax/σamax

plimsL = al -> [cbarmin(al), Tmin, 0.1, 0.1, 0.1]
plimsU = [1000.0, 1000.0, 1000.0, 100.0, 100.0]

# # One-off ID or opt ---------

# param1, _, traj1, unactErr = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.1, plimsL(1.6), plimsU; Fext_pdep=true, test=false, testTrajReconstruction=false, print_level=1, max_iter=200)
# display(param1')

# traj2 = cu.fixTrajWithDynConst(m, opt, traj1, param1)
# # cu.optAffine(m, opt, traj1, param1, 1, R_WTS, 0.1, cbarmin(1.5); Fext_pdep=false, test=true, print_level=2)
# pl1 = plotTrajs(m, opt, trajt, [param0, param1, param1], [traj0, traj1, traj2])
# plot(pl1...)

# # The actuator data does not correspond to the kinematics in any way (esp. without params)
# # # 1. Try to find the best params *assuming* these are the correct inputs. ID mode
# # param1, paramObj, traj1 = cu.optAffine(m, opt, traj0, param0, 2, (zeros(4,4), 0, 1.0*ones(1,1)), 0.3, cbarmin; Fext_pdep=false, test=false, print_level=1)

# # # 2. Try to optimize
# param1, paramObj, traj1 = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.3, cbarmin(0.5); Fext_pdep=false, test=false, print_level=1)

# # mwings = collect(0.1:0.1:2)
# # plot(mwings, paramObj.([[param0[1:2];mwing] for mwing in mwings]))

# display([param0, param1])
# pls = plotParamImprovement(m, opt, trajt, [param0, param1], [traj0, traj1], paramObj)
# plot(pls...)

# Debug components ----------------

function debugComponentsPlot(traj, param; optal=nothing)
	if !isnothing(optal)
		param1, _, traj1, unactErr = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.1, plimsL(optal), plimsU; Fext_pdep=true, test=false, testTrajReconstruction=false, print_level=1, max_iter=200)
	else
		param1 = param
		traj1 = traj
	end
	yo, HMqT, HC, Hg, Hgact, HF = cu.paramAffine(m, opt, traj, param, R_WTS; Fext_pdep=true, debugComponents=true)
	pt0, T0 = cu.getpt(m, param)
	inertial = zeros(2,N)
	stiffdamp = similar(inertial)
	stiffdampa = similar(inertial)
	aero = similar(inertial)

	for k=1:N
		inertial[:,k] = (HMqT(yo(k), yo(k+1)) - HMqT(yo(k), yo(k)) + opt.fixedδt * HC(yo(k))) * pt0
		stiffdamp[:,k] = (opt.fixedδt * Hg(yo(k))) * pt0
		stiffdampa[:,k] = (opt.fixedδt * Hgact(yo(k))) * pt0
		aero[:,k] = (opt.fixedδt * HF(yo(k))) * pt0
	end

	function plotComponents(i, ylbl)
		pl = plot(inertial[i,:], linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, stiffdamp[i,:], linewidth=2, label="g")
		plot!(pl, stiffdampa[i,:], linewidth=2, label="ga")
		plot!(pl, aero[i,:], linewidth=2, label="a")
		tot = inertial[i,:]+stiffdamp[i,:]+stiffdampa[i,:]+aero[i,:]
		plot!(pl, tot, linewidth=2, linestyle=:dash, label="tot")

		pl2 = plot(aero[i,:], linewidth=2, label="drag", legend=:outertopright)
		plot!(pl2, tot, linewidth=2, label="act")
		return pl, pl2
	end

	pl1 = plotTrajs(m, opt, trajt, [param], [traj])
	pls, plcomp = plotComponents(1, "stroke")
	plh, _ = plotComponents(2, "hinge")

	# Note that gamma is here
	println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2, ", optal = ", optal)

	return pl1..., pls, plh, plcomp
end
pls = debugComponentsPlot(traj0, param0; optal=0.8)
plot(pls..., size=(800,600))
gui()

error("TEST")

# many sims (scale) --------------

function maxuForMinAvgLift(al)
	print("plimsL[1:2] = ", plimsL(al)[1:2])
	param1, _, traj1, unactErr = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.1, plimsL(al), plimsU; Fext_pdep=true, test=false, testTrajReconstruction=false, print_level=1, max_iter=200)
	kΨ, bΨ = param1[4:5]
	return [param1; norm(traj1[(N+1)*ny:end], Inf); norm(unactErr, Inf)]
end

minlifts = 0.1:0.2:2.0
llabels = [
	"chord",
	"T",
	"mwing",
	"hinge k",
	"hinge b"
]

res = hcat(maxuForMinAvgLift.(minlifts)...)'
np = length(param0)
p1 = plot(minlifts, res[:,1:np], xlabel="min avg lift [mN]", label=llabels, ylabel="design params", linewidth=2, legend=:topleft)
p2 = plot(minlifts, res[:,np+1], xlabel="min avg lift [mN]", ylabel="umin [mN]", linewidth=2, legend=false)
p3 = plot(minlifts, res[:,np+2], xlabel="min avg lift [mN]", ylabel="unact err", linewidth=2, label="err", legend=false)
plot(p1, p2, p3)

# ! pick one
# res = maxuForMinAvgLift(3)
# ptest = res[1:np]

# res0 = maxuForMinAvgLift(avgLift0)
# res0[np+2]

# # traj opt ------------------------------------

# # εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# # prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj)

# # plot(plot(prob.g), plot(prob.mult_g), size=(900,400))

# mo = nothing#cu.paramoptQPSetup(m, opt, traj0; scaling=false, verbose=false)
# Q = nothing # FIXME: compile error

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
# pl1 = plotTrajs(m, opt, trajt, params, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# # plot(pl1...)
# # gui()
