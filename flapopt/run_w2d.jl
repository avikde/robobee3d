
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using BenchmarkTools, MAT, Dierckx
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
cycleFreqLims = [0.3,0.01]#[0.165,0.165]#[0.4, 0.03] # [KHz]
dtlims = 1.0 ./ (N*cycleFreqLims)
POPTS = cu.ParamOptOpts(
	τinds=[2,5], 
	R=(zeros(4,4), reshape([10.0],1,1), 0.0*I), # middle one is mech pow
	plimsL = [0.1, 1.0, 0.1, 0.5, 0, 20.0, dtlims[1]],
	plimsU = [400.0, 3.5, 100.0, 20.0, 100.0, 500.0, dtlims[2]],
	εunact = 1.0, # 0.1 default. Do this for now to iterate faster
	uinfnorm = true,
	unactWeight = 1.0
)
includet("w2d_shift.jl")
# FUNCTIONS GO HERE -------------------------------------------------------------

const SCALING1_FNAME = "scaling1.zip"
function scaling1(m::Wing2DOFModel, opt, traj, param, xs, minlifts, τ21ratiolim; kwargs...)
	np = length(param)
	function scaling1single(x, minlift)
		r = opt1(m, traj, param, 1, minlift, τ21ratiolim; Φ=x, kwargs...)
		return [x; minlift; r["param"]; r["u∞"]; r["al"]; r["δact"]; mean(abs.(r["mechPow"])); r["FD∞"]; norm(r["unactErr"], Inf)]
	end
	results = [scaling1single(x, minlift) for minlift in minlifts, x in xs] # reversed
	resdict = Dict(
		"Phis" => xs, "minlifts" => minlifts, "results" => results
	)
	# ^ returns a 2D array result arrays
	matwrite(SCALING1_FNAME, resdict; compress=true)

	return resdict
end

function scaling1disp(resarg; useFDasFact=true, scatterOnly=false, xpl=nothing, ypl=nothing, s=nothing, Fnom=75, mactline=7000)
	np = length(param0)
	resdict = typeof(resarg) == String ? matread(resarg) : resarg
	mactRobobee = Fnom*σamax

	# Produced unstructured xi,yi,zi data
	xi = Float64[]
	ARi = Float64[]
	Awi = Float64[]
	FLi = Float64[]
	Phii = Float64[]
	mli = Float64[]
	Lwi = Float64[]
	macti = Float64[]  # times Robobee act
	powi = Float64[]
	freqi = Float64[]
	Ti = Float64[]
	for res in resdict["results"]
		Phi = deg2rad(res[1])
		param = res[2+1:2+np]
		stats = res[2+np+1:end]
		Lw = param[6]/sqrt(param[1])

		# Append to unstructured data
		append!(Phii, Phi)
		append!(mli, res[2])
		append!(Lwi, Lw)
		append!(Awi, param[6])
		append!(ARi, Lw/sqrt(param[1]))
		append!(xi, Phi*Lw)
		append!(FLi, stats[2])
		append!(macti, stats[3] * (useFDasFact ? stats[5] : stats[1])/mactRobobee)
		append!(powi, stats[4])
		append!(freqi, 1000/(N*param[np]))
		append!(Ti, param[2])
	end

	# Output the plots
	pl1 = scatter(xlabel="Phi", ylabel="Lw", legend=false)
	pl2 = scatter(xlabel="x", ylabel="FL", legend=false)
	scatter!(pl1, Phii, Lwi)
	scatter!(pl2, xi, FLi)

	retpl = [pl1, pl2]

	if !scatterOnly
		# Plot range changes depending on opt results (see scatter)
		if isnothing(xpl)
			xpl = [minimum(xi), maximum(xi)]
			ypl = [minimum(yi), maximum(yi)]
		end
		X = range(xpl[1], xpl[2], length=50)
		Y = range(ypl[1], ypl[2], length=50)
		if isnothing(s)
			s = length(xi)
		end

		function contourFromUnstructured(xi, yi, zi; title="")
			# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
			# println("Total points = ", length(xi))
			spl = Spline2D(xi, yi, zi; s=s)
			ff(x,y) = spl(x,y)
			return contour(X, Y, ff, 
				titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
				xlabel="x [mm]", ylabel="FL [mg]", title=title,
				xlims=xpl, ylims=ypl)
		end

		# pl1 = plot(xs, [res[6]/res[1] for res in results], xlabel="Phi", ylabel="Lw", lw=2)

		plmact = contourFromUnstructured(xi, FLi, macti; title="mact [x Robobee]")
		plot!(plmact, X, mactline./X, lw=2, color=:black, ls=:dash, label="")

		append!(retpl, [plmact, 
			# scatter3d(xi, FLi, macti, camera=(90,40)),
			contourFromUnstructured(xi, FLi, powi; title="Avg mech pow [mW]"),
			contourFromUnstructured(xi, FLi, Awi; title="Aw [mm^2]"),
			contourFromUnstructured(xi, FLi, ARi; title="ARi"),
			# scatter3d(xi, FLi, powi, camera=(10,40)),
			# contourFromUnstructured(xi, FLi, rad2deg.(Phii); title="Phi"), 
			# contourFromUnstructured(xi, FLi, mli; title="ml"), 
			contourFromUnstructured(xi, FLi, freqi; title="freq [Hz]"), 
			contourFromUnstructured(xi, FLi, Ti; title="T1 [rad/mm]")])
	end
	
	return retpl
end

"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(ret, minlifts, τ21ratiolim; kwargs...)
	traj, param = ret["traj"], ret["param"]
	function maxuForMinAvgLift(al)
		r = opt1(m, traj, param, 1, al, τ21ratiolim; kwargs...)
		# kΨ, bΨ = param2[4:5]
		uu = r["traj"][(N+1)*ny:end]
		return [r["param"]; norm(uu, Inf); norm(r["unactErr"], Inf); norm(uu, 2)/N; r["al"]]
	end
	llabels = [
		"chord",
		"T1",
		"mwing",
		"hinge k",
		"hinge b",
		"T2",
		"Aw",
		"dt"
	]
	minliftsmg = minlifts

	res = hcat(maxuForMinAvgLift.(minlifts)...)'
	display(res)
	np = length(param0)
	actualliftsmg = res[:,np+4]
	p1 = plot(actualliftsmg, res[:,POPTS.τinds], xlabel="avg lift [mg]", label=llabels[POPTS.τinds], ylabel="T1,T2", linewidth=2, legend=:topleft)
	p2 = plot(actualliftsmg, [res[:,np+1]  res[:,np+3]], xlabel="avg lift [mg]", ylabel="umin [mN]", linewidth=2, legend=:topleft, label=["inf","2","al"])
	hline!(p2, [75], linestyle=:dash, color=:black, label="robobee act")
	p3 = plot(actualliftsmg, 1000 ./ (N*res[:,np]), xlabel="avg lift [mg]", ylabel="Cycle freq [Hz]", linewidth=2, legend=false)
	p4 = plot(actualliftsmg, res[:,np+2], xlabel="avg lift [mg]", ylabel="unact err", linewidth=2, legend=false)

	return p1, p2, p3, p4
end

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(ret, Tratios, minals)
	function maxu(τ21ratiolim, minal)
		rr = opt1(m, ret["traj"], ret["param"], 1, minal, τ21ratiolim)
		return [rr["u∞"]; rr["al"]; rr["FD∞"]; rr["param"]]
	end

	res = [[T;al;maxu(T,al)] for T in Tratios, al in minals]
	matwrite(NLBENEFIT_FNAME, Dict("res" => res); compress=true)
	return res
end

function plotNonlinBenefit(fname; s=100)
	results = matread(fname)["res"]
	
	# Row 1 has Tratio=0 (res[1,:])
	for r=size(results,1):-1:1
		for c=1:size(results,2)
			resT0 = results[1,c]
			# normalize
			results[r,c][3] /= resT0[3]
			results[r,c][5] /= resT0[5]
		end
	end
	xyzi = zeros(length(results[1,1]),0)
	for res in results
		xyzi = hcat(xyzi, res)
	end
	# lift to mg
	params = xyzi[6:end,:]

	xpl = [0,3]
	ypl = [140, 170]
	X = range(xpl[1], xpl[2], length=50)
	Y = range(ypl[1], ypl[2], length=50)

	function contourFromUnstructured(xi, yi, zi; title="")
		# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
		# println("Total points = ", length(xi))
		spl = Spline2D(xi, yi, zi; s=s)
		ff(x,y) = spl(x,y)
		return contour(X, Y, ff, 
			titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
			xlabel="T ratio", ylabel="FL [mg]", title=title,
			xlims=xpl, ylims=ypl)
	end
    
	return [
		# scatter(xyzi[1,:], xyzi[4,:]),
		# scatter3d(xyzi[1,:], xyzi[4,:], xyzi[3,:]),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], xyzi[3,:]; title="Nonlinear transmission benefit []"),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], params[2,:]; title="T1 [rad/mm]"),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], params[6,:]; title="Aw [mm^2]"),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
	]
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

# resdict = scaling1(m, opt, traj0, param0, collect(60.0:10.0:120.0), collect(150:20:350), 2) # SLOW
# pls = scaling1disp("scaling1_0.1e3_new.zip"; scatterOnly=false, xpl=[19,37], ypl=[160,350], s=500, useFDasFact=true, Fnom=50) # Found this by setting useFDasFact=false, and checking magnitudes
# plot(pls..., size=(1000,600), window_title="Scaling1", dpi=200)
# savefig("scaling1.png")
# gui()
# error("i")

# ID
ret1 = KINTYPE==1 ? Dict("traj"=>traj0, "param"=>param0) : opt1(m, traj0, param0, 2, 0.1, 0.0) # In ID force tau2=0

# 2. Try to optimize
ret2 = opt1(m, ret1["traj"], ret1["param"], 1, 150; Φ=120, Rpow=10)#; print_level=3, max_iter=10000)

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
