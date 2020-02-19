
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
	kbo = [30, 0],
	ma = 6,
	ka = 240,
	Amp = deg2rad.([90, 140]))
ny, nu = cu.dims(m)

function getInitialParams()
	# robobee scale
	return 75, [3.2^2,  # cbar2[mm^2] (area/R)^2
		2.6666, # τ1 (from 3333 rad/m, [Jafferis (2016)])
		0.55, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
		2.5, # wΨ [mm]
		0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
		54.4, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
		0.0758 # dt
	]
	# # for 120deg stroke
	# return 75, [3.5^2,  # cbar2[mm^2] (area/R)^2
	# 	3.3, # τ1 (from 3333 rad/m, [Jafferis (2016)])
	# 	0.7, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
	# 	1.8, # wΨ [mm]
	# 	0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
	# 	55, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
	# 	0.0758 # dt
	# ]
end
uampl, param0 = getInitialParams()
σamax = 0.3 # [mm] constant? for robobee actuators

include("w2d_paramopt.jl")

# IMPORTANT - load which traj here!!!
KINTYPE = 1
N, trajt, traj0, opt, Φ0 = initTraj(m, param0, KINTYPE; uampl=uampl)
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
	ΔySpikyBound = 0.05,
	pdes = zeros(7),
	pdesQ = [0.,0.,0.,0.,0.,0.,3e4],
	centralDiff = true
)

# ret1 = KINTYPE==1 ? Dict("traj"=>traj0, "param"=>param0) : opt1(m, traj0, param0, 2, 0.1, 0.0) # In ID force tau2=0
# Try to "fix" the initial traj by getting uinf so that it is more feasible
ret1 = opt1(m, traj0, param0, 1, 100)

# These are the actual lims
cycleFreqLims = [0.3,0.05] # [KHz] -- in order to avoid first order integration errors try to keep a small dt
dtlims = 1.0 ./ (N*cycleFreqLims)
POPTS.plimsL .= [0.1, 1.0, 0.1, 0.5, 0, 20.0, dtlims[1]]
POPTS.plimsU .= [50.0, 3.5, 100.0, 20.0, 100.0, 500.0, dtlims[2]]

# includet("w2d_shift.jl")
# includet("w2d_debug.jl")
# SCRIPT RUN STUFF HERE -----------------------------------------------------------------------

## Scaling1 -------

# includet("w2d_scaling1.jl")
# # resdict = scaling1(m, opt, traj0, param0, 2, range(60, 120, length=7), range(160, 420, length=7), range(1e3, 1e5, length=4)) # SLOW
# pls = scaling1disp("scaling1_3d4.zip"; scatterOnly=false, xpl=[27,33], ypl=[220,420], s=5e5, useFDasFact=true, Fnom=75, mactline=10e3) # Found this by setting useFDasFact=false, and checking magnitudes
# # plot(pls..., size=(1000,600), window_title="Scaling1")#, dpi=200)
# plot(pls[3], pls[4], pls[6], pls[7], size=(600,500))
# # savefig("scaling1.png")
# gui()

## 

# scaling2(m, opt, traj0, param0, 90, range(1e3, stop=1e5, length=5), 180, 2) # SLOW

## Try to optimize ---------

# includet("w2d_paramopt.jl")
# ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 250; Φ=90, Qdt=3e4)#, print_level=3)

## Bigbee ---------

# kact = 4x?. for bigbee it seems like to get down to that frequency, it wants huge wings
# "high power": Opt Φ=120, Qdt=50000.0, minal=400, τ2/1 lim=2.0 => 0, [17.199 3.423 0.893 6.224 1.899 68.796 0.065], fHz=191.5, al[mg]=386.5, u∞=252.8, FD∞=149.0, pow=37.1, J=320.5, AR=4.0, x=34.7
# ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 400; Φ=120, Qdt=5e4)

# "low power": Opt Φ=120, Qdt=0.0, minal=300, τ2/1 lim=2.0 => 0, [22.083 3.282 1.147 5.273 6.566 88.333 0.097], fHz=129.1, al[mg]=308.5, u∞=185.0, FD∞=146.2, pow=22.8, J=248.0, AR=4.0, x=39.4
ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 300; Φ=120, Qdt=1e3, Rpow=1e0)

## Param space convexity plot -----------------
# includet("w2d_pplots.jl")
# debugConvexity(m, opt, POPTS, ret1)

## DEBUG ----------

# includet("w2d_debug.jl")
# pls = debugΔYEffect(N, ny, ret2, ret1)
# # plot(pls..., size=(1000,600))
# plot(pls[end], size=(400,200))
# gui()

## Components ---------

# pls = debugComponentsPlot(m, opt, POPTS, ret2)
# plot(pls..., size=(800,600))
# gui()

## NONLIN BENEFIT -----------------

# includet("w2d_nlbenefit.jl")

# # openLoopPlotFinal(m, opt, param0)

# # results = nonlinBenefit("nonlinbig5.zip", ret1, range(0, 3, length=10), range(180, 360, length=6); τ2eq=true, tol=5e-2) # SLOW
# pls = plotNonlinBenefit("nonlinbig5.zip", [180,360]; s=100)
# plot(pls...)
# gui()

##
## ----------------
# pls = scaleParamsForlift(ret1, 0.6:0.2:2.0, 2)
# plot(pls...)
