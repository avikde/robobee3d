
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
	0.65#= 1.5 =#, #k output
	0, #b output
	6, # ma
	0, # ba
	250#= 0 =#, # ka
	true) # bCoriolis
ny, nu = cu.dims(m)
param0 = [5.411,  # cbar[mm] (area/R)
	18.681, # τ1 (from 3333 rad/m, R=17, [Jafferis (2016)])
	0.866, # mwing[mg]
	22.633, # kΨ [mN-mm/rad]
	12.037, # bΨ [mN-mm/(rad/ms)]
	0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
	0.109 # dt
]
σamax = 0.3 # [mm] constant? for robobee actuators

includet("w2d_paramopt.jl")

# IMPORTANT - load which traj here!!!
KINTYPE = 1
N, trajt, traj0, opt, avgLift0 = initTraj(KINTYPE)

# Param opt init
cycleFreqLims = [0.4, 0.03] # [KHz]
dtlims = 1.0 ./ (N*cycleFreqLims)
POPTS = cu.ParamOptOpts(
	τinds=[2,6], 
	R=(zeros(4,4), 0, 1.0*I), 
	plimsL = [0.1, 10, 0.1, 0.1, 0.1, 0, dtlims[1]],
	plimsU = [1000.0, 1000.0, 1000.0, 100.0, 100.0, 100.0, dtlims[2]],
	εunact = 1.0 # 0.1 default. Do this for now to iterate faster
)
includet("w2d_shift.jl")
# FUNCTIONS GO HERE -------------------------------------------------------------

"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(ret, minlifts, τ21ratiolim)
	traj, param = ret["traj"], ret["param"]
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
	p1 = plot(minliftsmg, res[:,POPTS.τinds], xlabel="min avg lift [mg]", label=llabels[POPTS.τinds], ylabel="T1,T2", linewidth=2, legend=:topleft)
	p2 = plot(minliftsmg, [res[:,np+1]  res[:,np+3]], xlabel="min avg lift [mg]", ylabel="umin [mN]", linewidth=2, legend=:topleft, label=["inf","2"])
	p3 = plot(minliftsmg, 1000 ./ (N*res[:,np]), xlabel="min avg lift [mg]", ylabel="Cycle freq [Hz]", linewidth=2, legend=false)
	p4 = plot(minliftsmg, res[:,np+2], xlabel="min avg lift [mg]", ylabel="unact err", linewidth=2, legend=false)

	return p1, p2, p3, p4
end

function plotNonlinBenefit(ret)
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

# SCRIPT RUN STUFF HERE -----------------------------------------------------------------------

# ID
ret1 = KINTYPE==1 ? Dict("traj"=>traj0, "param"=>param0) : opt1(traj0, param0, 2, 0.1, 0.0) # In ID force tau2=0

# 2. Try to optimize
ret2 = opt1(ret1["traj"], ret1["param"], 1, 1.6)#; print_level=3, max_iter=10000)
# ret3 = opt1(ret1["traj"], ret1["param"], 1, 1.0; print_level=3, max_iter=10000)
# traj3, param3, paramObj, _ = opt1(traj2, param2, 1, 1.3)
# paramObj2(p) = paramObj([p; zeros((N+1)*ny)])
# pls = plotParamImprovement(m, opt, [param1, param2, param3], [traj1, traj2, traj3], paramObj2)
# plot(pls...)

# testManyShifts(ret1, [0], 0.6)

# retTest = Dict("traj"=>ret2["traj"], "param"=>ret2["param"])
# retTest["param"][2]

# pl1 = plotTrajs(m, opt, listOfParamTraj(ret1, ret2)...)
# plot(pl1...)

# ---------
pls = debugComponentsPlot(ret2)
plot(pls..., size=(800,600))

# -----------------
# pls = plotNonlinBenefit() # SLOW
# plot(pls...)

# # ----------------
# pls = scaleParamsForlift(ret1, 0.6:0.2:1.6, 2)
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
# pl1 = plotTrajs(m, opt, params, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# # plot(pl1...)
# # gui()
