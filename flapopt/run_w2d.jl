
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using MAT, Dierckx
using Revise # while developing
import controlutils
cu = controlutils
include("w2d_model.jl")

# create an instance
# From Patrick ko = 300 mN-mm/rad but in input frame. This needs tuning...
# Noah said lumped stiffness is ~60% transmission and 40% actuator => k = 1.5 = ko + ka/T^2.
# For T=20, get ka = 240.
# To get ma, use the fact that actuator resonance is ~1KHz => equivalent ma = 240/(2*pi)^2 ~= 6mg
m = Wing2DOFModel(
	ko = 30.0,
	ma = 6,
	ka = 240,
	Amp = deg2rad.([90, 140]))
ny, nu = cu.dims(m)

function getInitialParams()
	# robobee scale
	return 75, [3.2^2,  # cbar2[mm^2] (area/R)^2
		2.6666, # τ1 (from 3333 rad/m, [Jafferis (2016)])
		0.73, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
		2.5, # wΨ [mm]
		0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
		54.4, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
		0.0758 # dt
	]
end
uampl, param0 = getInitialParams()
σamax = 0.3 # [mm] constant? for robobee actuators

include("w2d_paramopt.jl")

# IMPORTANT - load which traj here!!!
KINTYPE = 1
N, trajt, traj0, opt, Φ0 = initTraj(m, param0, KINTYPE; uampl=uampl)
# openLoopPlot(m, opt, param0; save=true)
avgLift0 = avgLift(m, opt, traj0, param0) # for minlift constraint
println("Avg lift initial [mg]=", round(avgLift0, digits=1))

# Param opt init
POPTS = cu.ParamOptOpts(
	τinds=[2,5], 
	plimsL = copy(param0),
	plimsU = copy(param0),
	R = (0.0*I, reshape([1e-2],1,1), 0.0*I, 1e4, 0.75), # Ryy, Ryu (mech pow), Ruu, wΔy, wlse
	εunact = 1.0, # 0.1 default. Do this for now to iterate faster
	εpoly = 1e-3,
	ΔySpikyBound = 0.03,
	pdes = zeros(7),
	pdesQ = [0.,0.,0.,0.,0.,0.,3e4],
	centralDiff = true
)

# ret1 = KINTYPE==1 ? Dict("traj"=>traj0, "param"=>param0) : opt1(m, traj0, param0, 2, 0.1, 0.0) # In ID force tau2=0
# Try to "fix" the initial traj by getting uinf so that it is more feasible
ret1 = opt1(m, traj0, param0, 1, 100)

# These are the actual lims
cycleFreqLims = [0.3,0.1] # [KHz] -- in order to avoid first order integration errors try to keep a small dt
dtlims = 1.0 ./ (N*cycleFreqLims)
POPTS.plimsL .= [0.1, 1.0, 0.1, 0.5, 0, 20.0, dtlims[1]]
POPTS.plimsU .= [50.0, 3.5, 100.0, 20.0, 100.0, 500.0, dtlims[2]]

# includet("w2d_shift.jl")
# includet("w2d_scaling1.jl")
# includet("w2d_nlbenefit.jl")
# includet("w2d_debug.jl")
# SCRIPT RUN STUFF HERE -----------------------------------------------------------------------

# # resdict = scaling1(m, opt, traj0, param0, collect(60.0:10.0:120.0), collect(150:20:350), 2) # SLOW
# pls = scaling1disp("scaling1_138.zip"; scatterOnly=false, xpl=[25,30], ypl=[200,400], s=500, useFDasFact=true, Fnom=50) # Found this by setting useFDasFact=false, and checking magnitudes
# plot(pls..., size=(1000,600), window_title="Scaling1", dpi=200)
# savefig("scaling1.png")
# gui()
# error("i")

# debug4()

# 2. Try to optimize
ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 180)#, print_level=3)
# pls = debugDeltaYEffect(N, ny, ret2)
# plot(pls..., size=(1000,600))

# ---------
pls = debugComponentsPlot(m, opt, POPTS, ret2)
plot(pls..., size=(800,600))

# # -----------------
# # nonlinBenefit(ret1, 0:0.5:3.0, 150:20:250) # SLOW
# pls = plotNonlinBenefit(NLBENEFIT_FNAME, [170,220]; s=20)
# plot(pls..., dpi=200)
# gui()
# savefig("nonlinbenefit_138.png")

# # ----------------
# pls = scaleParamsForlift(ret1, 0.6:0.2:2.0, 2)
# plot(pls...)
