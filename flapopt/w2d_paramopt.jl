
import controlutils
cu = controlutils
using Plots; gr()
include("Wing2DOF.jl")
include("LoadWingKinData.jl")

"""Produce initial traj
fix -- Make traj satisfy dyn constraint with these params?
"""
function initTraj(sim=false; fix=false, makeplot=false, Ψshift=0)
	if sim
		opt = cu.OptOptions(true, false, 0.2, 1, :none, 1e-8, false) # sim
		N = opt.boundaryConstraint == :symmetric ? 17 : 33
		trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0, 187)
	else
		# Load data
		opt = cu.OptOptions(true, false, 0.135, 1, :none, 1e-8, false) # real
		# N, trajt, traj0, lift, drag = loadAlignedData("data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.2445; strokeMult=m.R/(2*param0[2]), ForcePerVolt=0.8)
		N, trajt, traj0 = loadAlignedData("data/Bee1_Static_165Hz_180V_10KSF.mat", "data/Bee1_Static_165Hz_180V_7500sf.csv", 1250; sigi=1, strokeSign=1, strokeMult=m.R/(2*param0[2]), ForcePerVolt=75/100, vidSF=7320, Ψshift=Ψshift) # 75mN unidirectional at 200Vpp (from Noah)
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

"""One-off ID or opt"""
function opt1(traj, param, mode, minal, τ21ratiolim=2.0; testAffine=false, testAfter=false, testReconstruction=false, max_iter=4000, print_level=1)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	print("minal = ", minal, ", τ21ratiolim = ", τ21ratiolim, " => ")
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
	println(ret["param"]')
	return ret
end

listOfParamTraj(retlist) = [ret["param"] for ret in retlist], [ret["traj"] for ret in retlist]

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

	for k=1:N
		σo = yo(k)[1]
		inertial[:,k] = cu.Hτ(HMnc(yo(k), yo(k+1)) - HMnc(yo(k), yo(k)), σo) * pt0
		inertialc[:,k] = cu.Hτ(HMc(yo(k), yo(k+1)) - HMc(yo(k), yo(k)) + δt * HC(yo(k)), σo) * pt0
		stiffdamp[:,k] = cu.Hτ(δt * Hg(yo(k)), σo) * pt0
		stiffdampa[:,k] = cu.Hτ(δt * Hgact(yo(k)), σo) * pt0
		coriolis[:,k] = cu.Hτ(δt * HC(yo(k)), σo) * pt0
		aero[:,k] = cu.Hτ(δt * HF(yo(k)), σo) * pt0
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
		plot!(pl, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl, t2, tot / δt, linewidth=2, linestyle=:dash, label="tot")

		pl2 = plot(t2, aero[i,:] / δt, linewidth=2, label="-dr", legend=:outertopright)
		plot!(pl2, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl2, t2, coriolis[i,:] / δt, linewidth=2, label="cor")
		plot!(pl2, t2, inertiastiff / δt, linewidth=2, label="is")
		
		pl3 = plot(t2, inertial[i,:], linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, t2, inertialc[i,:], linewidth=2, label="ic")
		plot!(pl3, t2, inertial[i,:] + inertialc[i,:], linewidth=2, linestyle=:dash, label="itot")

		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, trajt, listOfParamTraj([ret])...)
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")

	# Note that gamma is here
	println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2)
	return pl1[[1,2,4,5]]..., pls, plh, plcomp, plis, plih
end
