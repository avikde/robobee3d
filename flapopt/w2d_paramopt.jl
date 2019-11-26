
import controlutils
cu = controlutils
using Plots; gr()
include("Wing2DOF.jl")
include("LoadWingKinData.jl")

"""From Chen (2016) IROS, find a usable estimate of the wing density to use as a bound."""
function estimateWingDensity(test=false)
	# Copy-pasted data from [Chen (2016)]. Each of these should presumably be a different "density" (instantiated by dspar).
	# For constraining the design, could use the lowest one of these.
	dspar = [0.14 0.17 0.20 0.23 0.26 0.29]'
	Ixx = [1.91 2.25 2.56 2.90 3.27 3.64]' # spanwise moment of inertia, i.e. in the stroke dynamics (through R^2)
	Izz = [40.6 48.8 57.2 65.6 73.4 82.7]' # chordwise moment of inertia, i.e. inertia in the hinge dynamics
	cbar0 = 5.345 # Use dC from LoadWingKinData.jl
	# Want to fit a density such that Ixx
	mwings = Ixx/(m.R/2)^2 #<- this TODO: this does not make much sense... check with Noah
	mwings2 = Izz/(2*cbar0^2*γ^2)
	meanpred = mwings2 ./ cbar0
	# FIXME: for robobee design this is producing torques that are too high
	# rholims = [meanpred[1] * cbar0, meanpred[end] * cbar0]
	rholims = [0.16, 0.35] # from param0 = 0.16

	if test
		p1 = plot(dspar, mwings, ylabel="mwing from Ixx", lw=2)
		p2 = plot(dspar, mwings2, ylabel="mwing from Izz", lw=2)
		hline!(p2, rholims)
		plot(p1, p2)
		gui()
		error("rholims", rholims, " from param0 ", param0[3]/param0[1])
	end
	
	# return param0[3]/param0[1] <- 0.16
	return rholims
end


"""Produce initial traj
kinType -- 0 => ID'ed real data, 1 => openloop sim with param0 then truncate, 2 => generate kinematics(t)
fix -- Make traj satisfy dyn constraint with these params?
"""
function initTraj(kinType=0; fix=false, makeplot=false, Ψshift=0)
	if kinType==1
		opt = cu.OptOptions(true, false, 0.135, 1, :none, 1e-8, false) # sim
		N = opt.boundaryConstraint == :symmetric ? 23 : 45
		trajt, traj0 = createInitialTraj(m, opt, N, 0.165, [1e3, 1e2], param0, 170; uampl=100)
	elseif kinType==0
		# Load data
		opt = cu.OptOptions(true, false, 0.135, 1, :none, 1e-8, false) # real
		# N, trajt, traj0, lift, drag = loadAlignedData("data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.2445; strokeMult=m.R/(2*param0[2]), ForcePerVolt=0.8)
		N, trajt, traj0 = loadAlignedData("data/Bee1_Static_165Hz_180V_10KSF.mat", "data/Bee1_Static_165Hz_180V_7500sf.csv", 1250; sigi=1, strokeSign=1, strokeMult=m.R/(2*param0[2]), ForcePerVolt=75/100, vidSF=7320, Ψshift=Ψshift) # 75mN unidirectional at 200Vpp (from Noah)
	else
		error("Not implemented")
	end

	if fix
		traj0 = cu.fixTrajWithDynConst(m, opt, traj0, param0)
	end

	if makeplot
		pl1 = plotTrajs(m, opt, [param0], [traj0])
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

"""One-off ID or opt"""
function opt1(traj, param, mode, minal, τ21ratiolim=2.0; testAffine=false, testAfter=false, testReconstruction=false, max_iter=4000, print_level=1)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	print(mode==2 ? "ID" : "Opt", " minal=", minal, ", τ2/1 lim=", τ21ratiolim, " => ")

	# Poly constraint
	rholims = estimateWingDensity()
	Cp = Float64[-1  0  0  0  0  0  0;
		0  -τ21ratiolim  0  0  0  1  0;
		rholims[1]   0  -1  0  0  0  0; # wing density mw >= cbar*ρ1
		-rholims[2]   0  1  0  0  0  0] # wing density mw <= cbar*ρ2
	cbarmin = minAvgLift -> param0[1] * minAvgLift / avgLift0
	dp = [-cbarmin(minal); 0; 0; 0]

	ret = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, Cp=Cp, dp=dp, print_level=print_level, max_iter=max_iter, testTrajReconstruction=testReconstruction)
	# append unactErr
	ret["unactErr"] = ret["eval_g"](ret["x"])[1:N] # 1 unact DOF
	ret["al"] = avgLift(m, opt, ret["traj"], ret["param"])
	uu = ret["traj"][(N+1)*ny:end]
	ret["u∞"] = norm(uu, Inf)
	if testAfter
		cu.affineTest(m, opt, ret["traj"], ret["param"], POPTS)
	end
	println(round.(ret["param"]', digits=3), ", fHz=", round(1000/(N*ret["param"][end]), digits=1), ", al[mg]=", round(ret["al"] * 1000/9.81, digits=1), ", u∞=", round(ret["u∞"], digits=1), ", J=", round(ret["eval_f"](ret["x"]), digits=1))
	return ret
end

listOfParamTraj(rets...) = [ret["param"] for ret in rets], [ret["traj"] for ret in rets]

"""Debug components in a traj"""
function debugComponentsPlot(ret)
	traj1, param1 = ret["traj"], ret["param"]
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, ret["traj"])

	# Get the components
	yo, HMnc, HMc, HC, Hg, Hgact, HF = cu.paramAffine(m, opt, traj1, param1, POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, param1)
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
		inertial[:,k] = [cu.Hτ(HMnc(yo(k), yo(k+1)) - HMnc(yo(k), yo(k)), σo)  H0] * pt0
		inertialc[:,k] = [cu.Hτ(HMc(yo(k), yo(k+1)) - HMc(yo(k), yo(k)), σo)  cu.Hτ(HC(yo(k)), σo)] * pt0
		stiffdamp[:,k] = [H0  cu.Hτ(Hg(yo(k)), σo)] * pt0
		stiffdampa[:,k] = [H0  cu.Hτ(Hgact(yo(k)), σo)] * pt0
		coriolis[:,k] = [H0  cu.Hτ(HC(yo(k)), σo)] * pt0
		aero[:,k] = [H0  cu.Hτ(HF(yo(k)), σo)] * pt0
	end

	# # get the instantaneous transmission ratio at time k
	# Tvec = [cu.transmission(m, yo(k), param1; o2a=true)[2] for k=1:N]

	t2 = trajt[1:end-1]

	function plotComponents(i, ylbl)
		# grav+inertial -- ideally they would "cancel" at "resonance"
		inertiastiff = inertial[i,:]+inertialc[i,:]+stiffdamp[i,:]+stiffdampa[i,:]
		tot = inertiastiff+aero[i,:]

		pl = plot(t2, (inertial[i,:] + inertialc[i,:]) / δt, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, t2, stiffdamp[i,:] / δt, linewidth=2, label="g")
		plot!(pl, t2, stiffdampa[i,:] / δt, linewidth=2, label="ga")
		# plot!(pl, t2, aero[i,:], linewidth=2, label="a")
		# plot!(pl, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl, t2, tot / δt, linewidth=2, linestyle=:dash, label="act")

		pl2 = plot(t2, aero[i,:] / δt, linewidth=2, label="-dr", legend=:outertopright)
		plot!(pl2, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl2, t2, coriolis[i,:] / δt, linewidth=2, label="cor")
		plot!(pl2, t2, inertiastiff / δt, linewidth=2, label="is")
		
		pl3 = plot(t2, inertial[i,:] / δt, linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, t2, inertialc[i,:] / δt, linewidth=2, label="ic")
		plot!(pl3, t2, (inertial[i,:] + inertialc[i,:]) / δt, linewidth=2, linestyle=:dash, label="itot")
		plot!(pl3, t2, -(stiffdamp[i,:] + stiffdampa[i,:]) / δt, linewidth=2, label="-gtot")

		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, listOfParamTraj(ret)...)
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")

	# Note that gamma is here
	# println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2)
	return pl1[[1,2,4,5]]..., pls, plh, plcomp, plis, plih
end
