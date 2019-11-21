
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
	plimsU = [1000.0, 1000.0, 1000.0, 100.0, 100.0, 100.0]
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
function opt1(traj, param, mode, minal, τ21ratiolim=2.0; testAffine=false, testAfter=false, testReconstruction=false, max_iter=4000, print_level=1)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	Cp = Float64[-1  0  0  0  0  0;
		0  -τ21ratiolim  0  0  0  1]
	cbarmin = minAvgLift -> param0[1] * minAvgLift / avgLift0
	dp = [-cbarmin(minal); 0]
	ret = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, Cp=Cp, dp=dp, print_level=print_level, max_iter=max_iter, testTrajReconstruction=testReconstruction)
	# append unactErr
	ret["unactErr"] = ret["eval_g"](ret["x"])[1:N] # 1 unact DOF
	if testAfter
		cu.affineTest(m, opt, ret["traj"], ret["param"], POPTS)
	end
	println("minal = ", minal, ", τ21ratiolim = ", τ21ratiolim, " => ", ret["param"]')
	return ret
end

"""Debug components in a traj"""
function debugComponentsPlot(ret)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, ret["traj"])

	# Get the components
	yo, HMnc, HMc, HC, Hg, Hgact, HF = cu.paramAffine(m, opt, ret["traj"], ret["param"], POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, ret["param"])
	inertial = zeros(2,N)
	inertialc = similar(inertial)
	stiffdamp = similar(inertial)
	stiffdampa = similar(inertial)
	aero = similar(inertial)

	for k=1:N
		σo = yo(k)[1]
		inertial[:,k] = cu.Hτ(HMnc(yo(k), yo(k+1)) - HMnc(yo(k), yo(k)), σo) * pt0
		inertialc[:,k] = cu.Hτ(HMc(yo(k), yo(k+1)) - HMc(yo(k), yo(k)) + δt * HC(yo(k)), σo) * pt0
		stiffdamp[:,k] = cu.Hτ(δt * Hg(yo(k)), σo) * pt0
		stiffdampa[:,k] = cu.Hτ(δt * Hgact(yo(k)), σo) * pt0
		aero[:,k] = cu.Hτ(δt * HF(yo(k)), σo) * pt0
	end

	# # get the instantaneous transmission ratio at time k
	# Tvec = [cu.transmission(m, yo(k), param1; o2a=true)[2] for k=1:N]

	function plotComponents(i, ylbl)
		pl = plot(inertial[i,:] + inertialc[i,:], linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, stiffdamp[i,:], linewidth=2, label="g")
		plot!(pl, stiffdampa[i,:], linewidth=2, label="ga")
		plot!(pl, aero[i,:], linewidth=2, label="a")
		tot = inertial[i,:]+inertialc[i,:]+stiffdamp[i,:]+stiffdampa[i,:]+aero[i,:]
		plot!(pl, tot, linewidth=2, linestyle=:dash, label="tot")

		pl2 = plot(aero[i,:] / δt, linewidth=2, label="-dr(af)", legend=:outertopright)
		plot!(pl2, traj1[(N+1)*ny+1:end], linewidth=2, label="actf")
		
		pl3 = plot(inertial[i,:], linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, inertialc[i,:], linewidth=2, label="ic")
		plot!(pl3, inertial[i,:] + inertialc[i,:], linewidth=2, linestyle=:dash, label="itot")
		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, trajt, [ret["param"]], [ret["traj"]])
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")

	# Note that gamma is here
	println("param = ", ret["param"]', ", Iw = ", ret["param"][3] * (0.5 * ret["param"][1])^2)
	return pl1[[1,2,4,5]]..., pls, plh, plcomp, plis, plih
end

"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(traj, param, minlifts, τ21ratiolim)
	function maxuForMinAvgLift(al)
		r = opt1(traj, param, 1, al, τ21ratiolim)
		# kΨ, bΨ = param2[4:5]
		uu = r["traj"][(N+1)*ny:end]
		return [r["param"]; norm(uu, Inf); norm(r["unactErr"], Inf); norm(uu, 2)/N]
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
	p1 = plot(minliftsmg, res[:,POPTS.τinds], xlabel="min avg lift [mg]", label=llabels[POPTS.τinds], ylabel="design params", linewidth=2, legend=:topleft)
	p2 = plot(minliftsmg, [res[:,np+1]  res[:,np+3]], xlabel="min avg lift [mg]", ylabel="umin [mN]", linewidth=2, legend=:topleft, label=["inf","2"])
	p3 = plot(minliftsmg, res[:,np+2], xlabel="min avg lift [mg]", ylabel="unact err", linewidth=2, label="err", legend=false)

	return p1, p2, p3
end

function plotNonlinBenefit()
    # First plot the param landscape
    pranges = [
        0:0.15:3.0, # τ21ratiolim
        0.5:0.1:1.5 # minal
    ]
    labels = [
        "nonlin ratio",
        "minal"
	]
	
	function maxu(τ21ratiolim, minal)
		r = opt1(traj1, param1, 1, minal, τ21ratiolim)
		return norm(r["traj"][(N+1)*ny+1:end], Inf)
	end

	function plotSlice(i1, i2)
		zgrid = [maxu(x,y) for y in pranges[i2], x in pranges[i1]] # reversed: see https://github.com/jheinen/GR.jl/blob/master/src/GR.jl
		# to get the improvement, divide each metric by the performance at τ2=0
		maxuatτ2_0 = zgrid[:,1]
		zgrid = zgrid ./ repeat(maxuatτ2_0, 1, length(pranges[i1]))

		pl = contourf(pranges[i1], pranges[i2], zgrid, fill=true, seriescolor=cgrad(:bluesreds), xlabel=labels[i1], ylabel=labels[i2])
        # just in case
        xlims!(pl, (pranges[i1][1], pranges[i1][end]))
        ylims!(pl, (pranges[i2][1], pranges[i2][end]))
        return pl
    end
    
    return (plotSlice(1, 2),)
end

# Test feasibility
function paramTest(p, paramConstraint)
	xtest = [p; zeros((N+1)*ny)]
	g = paramConstraint(xtest)
	gunact = g[1:(N+1)]
	grest = g[(N+1)+1:end]
	# have [gpolycon (2); gtransmission (1); ginfnorm (2*N)]
	println("Feas: should be nonpos: ", maximum(grest), "; unact: ", maximum(abs.(gunact)) ,"; transmission: ", g[3])
	unew = cu.getTrajU(m, opt, traj1, p, POPTS)
	println("Obj: ", paramObj(xtest))
end

listOfParamTraj(retlist) = [ret["param"] for ret in retlist], [ret["traj"] for ret in retlist]

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
ret1 = opt1(traj0, param0, 2, 0.1)
# pl1 = plotTrajs(m, opt, trajt, [param1, param1], [traj0, traj1])
# plot(pl1...)

# 2. Try to optimize
ret2 = opt1(ret1["traj"], ret1["param"], 1, 0.8; print_level=3, max_iter=10000)
ret3 = opt1(ret1["traj"], ret1["param"], 1, 1.0; print_level=3, max_iter=10000)
# traj3, param3, paramObj, _ = opt1(traj2, param2, 1, 1.3)
# pl1 = plotTrajs(m, opt, trajt, [param1, param1, param2, param3], [traj0, traj1, traj2, traj3])
# # plot(pl1...)
# paramObj2(p) = paramObj([p; zeros((N+1)*ny)])
# pls = plotParamImprovement(m, opt, trajt, [param1, param2, param3], [traj1, traj2, traj3], paramObj2)
# plot(pls...)

pl1 = plotTrajs(m, opt, trajt, listOfParamTraj([ret2, ret3])...)
plot(pl1...)

# ---------
# pls = debugComponentsPlot(traj2, param2)
# plot(pls..., size=(800,600))

# -----------------
# pls = plotNonlinBenefit() # SLOW
# plot(pls...)

# ----------------
# pls = scaleParamsForlift(ret1["traj"], ret1["param"], 0.2:0.2:1.6, 0)
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
