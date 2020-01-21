
import controlutils
cu = controlutils
using ForwardDiff
using Plots; gr()
include("Wing2DOF.jl")
include("LoadWingKinData.jl")

# Indices into param array
const cb_idx = 1
const mw_idx = 3
const Aw_idx = 6

"""From Chen (2016) IROS, find a usable estimate of the wing density to use as a bound."""
function estimateWingDensity(test=false)
	# Copy-pasted data from [Chen (2016)]. Each of these should presumably be a different "density" (instantiated by dspar).
	# For constraining the design, could use the lowest one of these.
	dspar = [0.14 0.17 0.20 0.23 0.26 0.29]'
	Ixx = [1.91 2.25 2.56 2.90 3.27 3.64]' # spanwise moment of inertia, i.e. in the stroke dynamics (through R^2)
	Izz = [40.6 48.8 57.2 65.6 73.4 82.7]' # chordwise moment of inertia, i.e. inertia in the hinge dynamics
	cbar0 = 5.345 # Use dC from LoadWingKinData.jl
	# # Want to fit a density such that Ixx
	# mwings = Ixx/(m.R/2)^2 #<- this TODO: this does not make much sense... check with Noah
	# mwings2 = Izz/(2*cbar0^2*γ^2)
	# meanpred = mwings2 ./ cbar0
	# FIXME: for robobee design this is producing torques that are too high
	# rholims = [meanpred[1] * cbar0, meanpred[end] * cbar0]
	rholims = [0.01, 0.04] # from param0 = 0.52/54.4

	if test
		p1 = plot(dspar, mwings, ylabel="mwing from Ixx", lw=2)
		p2 = plot(dspar, mwings2, ylabel="mwing from Izz", lw=2)
		hline!(p2, rholims)
		plot(p1, p2)
		gui()
		error("rholims", rholims, " from param0 ", param0[mw_idx]/param0[Aw_idx])
	end
	
	# return param0[3]/param0[1] <- 0.16
	return rholims
end

"""Produce initial traj
kinType -- 0 => ID'ed real data, 1 => openloop sim with param0 then truncate, 2 => generate kinematics(t)
fix -- Make traj satisfy dyn constraint with these params?
"""
function initTraj(kinType=0; fix=false, makeplot=false, Ψshift=0, uampl=65, starti=172)
	if kinType==1
		opt = cu.OptOptions(false, false, 0.135, 1, :none, 1e-8, false) # sim
		N = opt.boundaryConstraint == :symmetric ? 23 : 45
		trajt, traj0 = createInitialTraj(m, opt, N, 0.165, [1e3, 1e2], param0, starti; uampl=uampl)
	elseif kinType==0
		# Load data
		cbar, τ1, mwing, wΨ, τ2, Aw, dt = param0
		R = Aw/cbar
		opt = cu.OptOptions(false, false, 0.135, 1, :none, 1e-8, false) # real
		# N, trajt, traj0, lift, drag = loadAlignedData("data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.2445; strokeMult=m.R/(2*param0[2]), ForcePerVolt=0.8)
		N, trajt, traj0 = loadAlignedData("data/Bee1_Static_165Hz_180V_10KSF.mat", "data/Bee1_Static_165Hz_180V_7500sf.csv", 1250; sigi=1, strokeSign=1, strokeMult=R/(2*τ1), ForcePerVolt=75/100, vidSF=7320, Ψshift=Ψshift) # 75mN unidirectional at 200Vpp (from Noah)
	else
		error("Not implemented")
	end

	if fix
		traj0 = cu.fixTrajWithDynConst(m, opt, traj0, param0)
	end

	if makeplot
		pl1 = plotTrajs(m, opt, [param0], [traj0]; legends=false)
		# pl1 = compareTrajToDAQ(m, opt, trajt, param0, traj0, lift, drag)
		plot(pl1...)
		gui()
		error("Initial traj")
	end
	
	# Constraint on cbar placed by minAvgLift
	avgLift0 = avgLift(m, opt, traj0, param0)
	println("Avg lift initial [mN]=", round(avgLift0, digits=4), ", in [mg]=", round(avgLift0*1000/9.81, digits=1))

	return N, trajt, traj0, opt, avgLift0
end

"""Linear approx of wing AR constraint at cbar https://github.com/avikde/robobee3d/issues/113. Returns a,b s.t. dot(a, [cb,Aw]) <= b is the constraint"""
function wingARconstraintLin(cbar; maxAR=4)
	# AR<=4: see https://github.com/avikde/robobee3d/pull/105#issuecomment-562761586 originally
	fwingAR(cbAw::Vector) = -maxAR*cbAw[1]^2+cbAw[2]
	ptang(cb) = [cb, -fwingAR([cb, 0])]
	p0 = ptang(cbar)
	Dfwing(x) = ForwardDiff.gradient(fwingAR, x)
	return Dfwing(p0), dot(Dfwing(p0),p0)
end

"""Linear approx of wing AR constraint at cbar https://github.com/avikde/robobee3d/issues/113. Returns a,b s.t. dot(a, [cb,Aw]) <= b is the constraint"""
function minLiftConstraintLin(minlift, param0, avgLift0)
	# Lift ~ Aw/dt^2 => fcons = Awdt2min - Aw/dt^2
	# Awdt2min = param0[Aw_idx]/param0[end]^2 * minlift / avgLift0
	Awdt2min = param0[Aw_idx] * minlift / avgLift0
	# println("hehe", Awdt2min)
	fminlift(Awdt::Vector) = -Awdt[1]#/Awdt[2]^2
	p0 = [param0[Aw_idx], param0[end]]
	Dfminlift(x) = ForwardDiff.gradient(fminlift, x)
	mla, mlb = Dfminlift(p0), -Awdt2min#dot(Dfminlift(p0),p0)
	# println("HI ", avgLift0, mla, mlb)
	return mla, mlb
end

"""One-off ID or opt"""
function opt1(traj, param, mode, minal, τ21ratiolim=2.0; testAffine=false, testAfter=false, testReconstruction=false, max_iter=4000, print_level=1, wARconstraintLinCbar=4.0)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	print(mode==2 ? "ID" : "Opt", " minal=", minal, ", τ2/1 lim=", τ21ratiolim, " => ")

    # cbar, τ1, mwing, wΨ, τ2, Aw, dt  = param
	# Poly constraint
	rholims = estimateWingDensity()
	wARa, wARb = wingARconstraintLin(wARconstraintLinCbar; maxAR=4)
	mla, mlb = minLiftConstraintLin(minal, param0, avgLift0)
	# Polytope constraint
	Cp = Float64[0  0  0  0  0  mla[1]  mla[2]; # min lift
		0  -τ21ratiolim  0  0  1  0  0; # transmission nonlinearity τ2 <= τ21ratiolim * τ1
		0   0  -1  0  0  rholims[1]  0; # wing density mw >= Aw*ρ1
		0   0  1  0  0  -rholims[2]  0; # wing density mw <= Aw*ρ2
		wARa[1]   0  0  0  0  wARa[2]  0] # wing AR
	dp = [mlb; 0; 0; 0; wARb]

	ret = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, Cp=Cp, dp=dp, print_level=print_level, max_iter=max_iter, testTrajReconstruction=testReconstruction)
	# append unactErr
	ret["unactErr"] = ret["eval_g"](ret["x"])[1:N] # 1 unact DOF
	ret["al"] = avgLift(m, opt, ret["traj"], ret["param"])
	uu = ret["traj"][(N+1)*ny:end]
	ret["u∞"] = norm(uu, Inf)
	if testAfter
		cu.affineTest(m, opt, ret["traj"], ret["param"], POPTS)
	end
	println(ret["status"], ", ", round.(ret["param"]', digits=3), ", fHz=", round(1000/(N*ret["param"][end]), digits=1), ", al[mg]=", round(ret["al"] * 1000/9.81, digits=1), ", u∞=", round(ret["u∞"], digits=1), ", J=", round(ret["eval_f"](ret["x"]), digits=1), ", AR=", round(ret["param"][Aw_idx]/ret["param"][cb_idx]^2, digits=1))
	return ret
end

listOfParamTraj(rets...) = [ret["param"] for ret in rets], [ret["traj"] for ret in rets]

"""Debug components in a traj"""
function debugComponentsPlot(m, opt, POPTS, ret)
	traj1, param1 = ret["traj"], ret["param"]
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, ret["traj"])

	# Get the components
	yo, HMnc, HMc, HC, Hg, Hgact, HF = cu.paramAffine(m, opt, traj1, param1, POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, param1)
	dt = param1[end] # need this for H after #110
	inertial = zeros(2,N)
	inertialc = similar(inertial)
	coriolis = similar(inertial)
	stiffdamp = similar(inertial)
	stiffdampa = similar(inertial)
	aero = similar(inertial)
	# to fill in for (non) inertial half of H https://github.com/avikde/robobee3d/pull/102
	H0 = zeros(size(inertial,1),length(pt0)÷2)

	for k=1:N
		σo = yo(k)[1]
		inertial[:,k] = [cu.Hτ(HMnc(yo(k), yo(k+1)) - HMnc(yo(k), yo(k)), σo)*dt  H0] * pt0
		inertialc[:,k] = [cu.Hτ(HMc(yo(k), yo(k+1)) - HMc(yo(k), yo(k)), σo)*dt  cu.Hτ(HC(yo(k)), σo)] * pt0
		stiffdamp[:,k] = [H0  cu.Hτ(Hg(yo(k)), σo)] * pt0
		stiffdampa[:,k] = [H0  cu.Hτ(Hgact(yo(k)), σo)] * pt0
		coriolis[:,k] = [H0  cu.Hτ(HC(yo(k)), σo)] * pt0
		aero[:,k] = [H0  cu.Hτ(HF(yo(k)), σo)] * pt0
	end

	# # get the instantaneous transmission ratio at time k
	# Tvec = [cu.transmission(m, yo(k), param1; o2a=true)[2] for k=1:N]

	t2 = collect(0:(N-1))*dt

	function plotComponents(i, ylbl)
		# grav+inertial -- ideally they would "cancel" at "resonance"
		inertiastiff = inertial[i,:]+inertialc[i,:]+stiffdamp[i,:]+stiffdampa[i,:]
		tot = inertiastiff+aero[i,:]

		pl = plot(t2, (inertial[i,:] + inertialc[i,:]) / dt, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, t2, stiffdamp[i,:] / dt, linewidth=2, label="g")
		plot!(pl, t2, stiffdampa[i,:] / dt, linewidth=2, label="ga")
		# plot!(pl, t2, aero[i,:], linewidth=2, label="a")
		plot!(pl, t2, traj1[(N+1)*ny+1:end], linewidth=2)
		plot!(pl, t2, tot/dt, linewidth=2, linestyle=:dash, label="act") # checked; this matches the row above ^

		pl2 = plot(t2, aero[i,:] / dt, linewidth=2, label="-dr", legend=:outertopright)
		plot!(pl2, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl2, t2, coriolis[i,:] / dt, linewidth=2, label="cor")
		plot!(pl2, t2, inertiastiff / dt, linewidth=2, label="is")
		
		pl3 = plot(t2, inertial[i,:] / dt, linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, t2, inertialc[i,:] / dt, linewidth=2, label="ic")
		plot!(pl3, t2, (inertial[i,:] + inertialc[i,:]) / dt, linewidth=2, linestyle=:dash, label="itot")
		plot!(pl3, t2, -(stiffdamp[i,:] + stiffdampa[i,:]) / dt, linewidth=2, label="-gtot")

		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, listOfParamTraj(ret)...)
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")

	# Note that gamma is here
	# println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2)
	return pl1[[1,2,4,6]]..., pls, plh, plcomp, plis, plih
end
