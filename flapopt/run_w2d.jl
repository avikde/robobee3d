
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using MAT, Dierckx
using Revise # while developing
import controlutils
cu = controlutils
include("Wing2DOF.jl")

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
	R = (0.0*I, reshape([1e-2],1,1), 0.0*I, 1e4, 1e-3, 1e0, 0), # Ryy, Ryu (mech pow), Ruu, wΔy, wu∞, wlse, wunact
	εunact = 1.0, # 0.1 default. Do this for now to iterate faster
	uinfnorm = false
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
# FUNCTIONS GO HERE -------------------------------------------------------------

# # Test feasibility
# function paramTest(p, paramConstraint)
# 	xtest = [p; zeros((N+1)*ny)]
# 	g = paramConstraint(xtest)
# 	gunact = g[1:(N+1)]
# 	grest = g[(N+1)+1:end]
# 	# have [gpolycon (2); gtransmission (1); ginfnorm (2*N)]
# 	println("Feas: should be nonpos: ", maximum(grest), "; unact: ", maximum(abs.(gunact)) ,"; transmission: ", g[3])
# 	unew = cu.getTrajU(m, opt, traj1, p, POPTS)
# 	println("Obj: ", paramObj(xtest))
# end

"See https://github.com/avikde/robobee3d/pull/136"
function debugDeltaYEffect(rr)
	pt, Hk, B, Js, actVec = rr["eval_f"](rr["x"]; debug=true)
	println("Js ", Js)
	dy = rr["x"][length(param0)+1:end]
	dely(k) = dy[(k-1)*ny+1:(k)*ny]
	dely0 = zeros(ny)
	unew = vcat([B' * Hk(k,dely(k),dely(k+1))[1] * pt for k=1:N]...)
	unew0 = vcat([B' * Hk(k,dely0,dely0)[1] * pt for k=1:N]...)
	p1 = plot([rr["traj"][(N+1)*ny+1:end]  actVec[:,1]  unew], lw=2, ls=[:solid :solid :dash])
	plot!(p1, unew0, lw=2, ls=:dash)
	return p1, plot(plot(dy[1:ny:(N+1)*ny]),
		plot(dy[2:ny:(N+1)*ny]),
		plot(dy[3:ny:(N+1)*ny]),
		plot(dy[4:ny:(N+1)*ny]))
end

"https://github.com/avikde/robobee3d/pull/137"
function debugGradient(rr)
	x = copy(rr["x"])
	x1s = 5.0:40
	fx1(x1, a=true) = rr["eval_f"]([x1; x[2:end]]; auto=a)
	function ddx1(x1, a=true)
		grad_f = similar(x)
		rr["eval_grad_f"]([x1; x[2:end]], grad_f; auto=a)
		return grad_f[1]
	end
	plot(
		plot([fx1.(x1s) fx1.(x1s, false)],lw=2,ls=[:solid :dash]), 
		plot([ddx1.(x1s)  ddx1.(x1s, false)],lw=2,ls=[:solid :dash])
	)
	gui()
	error("h")
end

# SCRIPT RUN STUFF HERE -----------------------------------------------------------------------

# # resdict = scaling1(m, opt, traj0, param0, collect(60.0:10.0:120.0), collect(150:20:350), 2) # SLOW
# pls = scaling1disp("scaling1_u2norm.zip"; scatterOnly=false, xpl=[30,35], ypl=[250,400], s=500, useFDasFact=true, Fnom=50) # Found this by setting useFDasFact=false, and checking magnitudes
# plot(pls..., size=(1000,600), window_title="Scaling1", dpi=200)
# savefig("scaling1.png")
# gui()
# error("i")

# ff = p -> cu.getpt(m, p)

# error("hi", ff(param0))

# 2. Try to optimize
ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 180; print_level=3#= , max_iter=10000 =#)
pls = debugDeltaYEffect(ret2)
plot(pls...)
gui()
error("hi")
# debugGradient(ret2)
# testManyShifts(ret1, [0], 0.6)

# retTest = Dict("traj"=>ret2["traj"], "param"=>ret2["param"])
# retTest["param"][2]

# pl1 = plotTrajs(m, opt, listOfParamTraj(ret1, ret2)...)
# plot(pl1...)

# ---------
pls = debugComponentsPlot(m, opt, POPTS, ret2)
plot(pls..., size=(800,600))

# # -----------------
# # nonlinBenefit(ret1, 0:0.3:3.0, 1.6:0.2:2.6) # SLOW
# pls = plotNonlinBenefit(NLBENEFIT_FNAME, s=500)
# plot(pls..., dpi=200)

# # ----------------
# pls = scaleParamsForlift(ret1, 0.6:0.2:2.0, 2)
# plot(pls...)
