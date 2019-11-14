
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
	0.55#= 1.5 =#, #k output
	0, #b output
	6, # ma
	0, # ba
	250#= 0 =#) # ka
ny, nu = cu.dims(m)
param0 = [3.2,  # cbar[mm] (area/R)
	28.33, # τ1 (from 3333 rad/m, R=17, [Jafferis (2016)])
	0.52, # mwing[mg]
	5, # kΨ [mN-mm/rad]
	3, # bΨ [mN-mm/(rad/ms)]
	0 # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
]

POPTS = cu.ParamOptOpts(
	τinds=[2,6], 
	R=(zeros(4,4), 0, 1.0*I), 
	plimsL = [0.1, 10, 0.1, 0.1, 0.1, 0],
	plimsU = [1000.0, 1000.0, 1000.0, 100.0, 100.0, 25.0]
)
σamax = 0.3 # [mm] constant? for robobee actuators

# FUNCTIONS GO HERE -------------------------------------------------------------

"""Produce initial traj
fix -- Make traj satisfy dyn constraint with these params?
"""
function initTraj(sim=false; fix=false, makeplot=false)
	if sim
		opt = cu.OptOptions(true, false, 0.2, 1, :none, 1e-8, false) # sim
		N = opt.boundaryConstraint == :symmetric ? 17 : 33
		trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0, 187)
	else
		# Load data
		opt = cu.OptOptions(true, false, 0.135, 1, :none, 1e-8, false) # real
		# N, trajt, traj0, lift, drag = loadAlignedData("data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.2445; strokeMult=m.R/(2*param0[2]), ForcePerVolt=0.8)
		N, trajt, traj0 = loadAlignedData("data/Bee1_Static_165Hz_180V_10KSF.mat", "data/Bee1_Static_165Hz_180V_7500sf.csv", 1250; sigi=1, strokeSign=1, strokeMult=m.R/(2*param0[2]), ForcePerVolt=75/100, vidSF=7320) # 75mN unidirectional at 200Vpp (from Noah)
	end

	if fix
		traj0 = cu.fixTrajWithDynConst(m, opt, traj0, param0)
	end

	if makeplot
		pl1 = plotTrajs(m, opt, trajt, [param0], [traj0])
		# pl1 = compareTrajToDAQ(m, opt, trajt, param0, traj0, lift, drag)
		plot(pl1...)
		gui()
	end

	return N, trajt, traj0, opt
end

# IMPORTANT - load which traj here!!!
N, trajt, traj0, opt = initTraj()

# Constraint on cbar placed by minAvgLift
avgLift0 = avgLift(m, opt, traj0, param0)

"""One-off ID or opt"""
function opt1(traj, param, mode, minal; testAffine=false, testAfter=false)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin
	Cp = Float64[-1  0  0  0  0  0]
	cbarmin = minAvgLift -> param0[1] * minAvgLift / avgLift0
	dp = [-cbarmin(minal)]
	param1, paramObj, traj1, unactErr = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, Cp=Cp, dp=dp, print_level=1, max_iter=4000)
	if testAfter
		cu.affineTest(m, opt, traj1, param1)
	end
	return traj1, param1, paramObj, unactErr
end

"""Debug components in a traj"""
function debugComponentsPlot(traj, param; optal=nothing)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
	optoptions = oaOpts(optal)

	if !isnothing(optal)
		traj1, param1, _, _ = opt1(traj, param, 1, optal)
	else
		param1 = param
		traj1 = traj
	end

	# Get the components
	yo, HMnc, HMc, HC, Hg, Hgact, HF = cu.paramAffine(m, opt, traj1, param1, POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, param1)
	inertial = zeros(2,N)
	inertialc = similar(inertial)
	stiffdamp = similar(inertial)
	stiffdampa = similar(inertial)
	aero = similar(inertial)

	for k=1:N
		inertial[:,k] = (HMnc(yo(k), yo(k+1)) - HMnc(yo(k), yo(k))) * pt0
		inertialc[:,k] = (HMc(yo(k), yo(k+1)) - HMc(yo(k), yo(k)) + δt * HC(yo(k))) * pt0
		stiffdamp[:,k] = (δt * Hg(yo(k))) * pt0
		stiffdampa[:,k] = (δt * Hgact(yo(k))) * pt0
		aero[:,k] = (δt * HF(yo(k))) * pt0
	end

	function plotComponents(i, ylbl)
		pl = plot(inertial[i,:] + inertialc[i,:], linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, stiffdamp[i,:], linewidth=2, label="g")
		plot!(pl, stiffdampa[i,:], linewidth=2, label="ga")
		plot!(pl, aero[i,:], linewidth=2, label="a")
		tot = inertial[i,:]+inertialc[i,:]+stiffdamp[i,:]+stiffdampa[i,:]+aero[i,:]
		plot!(pl, tot, linewidth=2, linestyle=:dash, label="tot")

		pl2 = plot(aero[i,:] * Tnew / δt, linewidth=2, label="-dr(af)", legend=:outertopright)
		plot!(pl2, traj1[(N+1)*4+1:end], linewidth=2, label="actf")
		
		pl3 = plot(inertial[i,:], linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, inertialc[i,:], linewidth=2, label="ic")
		plot!(pl3, inertial[i,:] + inertialc[i,:], linewidth=2, linestyle=:dash, label="itot")
		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, trajt, [param], [traj])
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")

	# Note that gamma is here
	println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2, ", optal = ", (!isnothing(optal) ? optal : "-"))
	return pl1[[1,2,4,5]]..., pls, plh, plcomp, plis, plih
end

"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(traj, param, minlifts)
	function maxuForMinAvgLift(al)
		println("plimsL[1:2] = ", plimsL(al)[1:2])
		traj1, param1, _, unactErr = opt1(traj, param, 1, al)
		kΨ, bΨ = param1[4:5]
		return [param1; norm(traj1[(N+1)*ny:end], Inf); norm(unactErr, Inf)]
	end
	llabels = [
		"chord",
		"T1",
		"mwing",
		"hinge k",
		"hinge b",
		"T2"
	]
	minliftsmg = minlifts .* 1000/9.81

	res = hcat(maxuForMinAvgLift.(minlifts)...)'
	display(res)
	np = length(param0)
	p1 = plot(minliftsmg, res[:,1:np], xlabel="min avg lift [mg]", label=llabels, ylabel="design params", linewidth=2, legend=:topleft)
	p2 = plot(minliftsmg, res[:,np+1], xlabel="min avg lift [mg]", ylabel="umin [mN]", linewidth=2, legend=false)
	p3 = plot(minliftsmg, res[:,np+2], xlabel="min avg lift [mg]", ylabel="unact err", linewidth=2, label="err", legend=false)

	return p1, p2, p3
end

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

# SCRIPT RUN STUFF HERE -----------------------------------------------------------------------

# ID
traj1, param1, paramObj, _ = opt1(traj0, param0, 2, 0.1)
# display(param1')
# pl1 = plotTrajs(m, opt, trajt, [param1, param1], [traj0, traj1])
# plot(pl1...)
# # 2. Try to optimize
# traj2, param2, paramObj, _ = opt1(traj1, param1, 1, 0.8)
# # display(param2')
# traj3, param3, paramObj, _ = opt1(traj2, param2, 1, 1.3
# # pl1 = plotTrajs(m, opt, trajt, [param1, param1, param2, param3], [traj0, traj1, traj2, traj3])
# # plot(pl1...)
# pls = plotParamImprovement(m, opt, trajt, [param1, param2, param3], [traj1, traj2, traj3], paramObj)
# plot(pls...)

# TEST manual params
Hk, yo, umeas, B, N = cu.paramAffine(m, opt, traj1, param1, POPTS, 1.0)
Δy0 = zeros((N+1)*ny)
testp(pnew) = cu.reconstructTrajFromΔy(m, opt, traj1, yo, Hk, B, Δy0, pnew)
pp = [8.44463,  13.9429,  0.235645,  23.9639,    8.29057,   0.0]
traj2 = testp(pp)
pptest = [8.44463,  11.0,  0.235645,  23.9639,    8.29057,   40.0]
traj3 = testp(pptest)
pl1 = plotTrajs(m, opt, trajt, [param1, param1, pp, pptest], [traj0, traj1, traj2, traj3])
plot(pl1...)

# ---------

# pls = debugComponentsPlot(traj1, param1; optal=1.2)
# plot(pls..., size=(800,600))
# gui()

# error("TEST")

# ----------------

# pls = scaleParamsForlift(traj1, param1, 0.5:0.2:1.5)
# plot(pls...)

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
# pl1 = plotTrajs(m, opt, trajt, params, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# # plot(pl1...)
# # gui()
