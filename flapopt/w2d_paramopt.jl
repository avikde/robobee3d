
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
	rholims = [0.013, 0.02] # from param0 = 0.7/54.4

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
function initTraj(m, param0, kinType=0; fix=false, makeplot=false, Ψshift=0, uampl=65, starti=214, verbose=true, freq=0.165, N=80)
	if kinType==1
		opt = cu.OptOptions(false, false, 1/(N*freq), 1, :none, 1e-8, false) # sim
		N = opt.boundaryConstraint == :symmetric ? N÷2 : N
		trajt, traj0 = createInitialTraj(m, opt, N, freq, [1e3, 1e2], param0, starti; uampl=uampl, verbose=verbose)
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

	@views tcoord(i) = traj0[i:ny:(N+1)*ny]
	currentAmpl(i) = maximum(tcoord(i)) - minimum(tcoord(i))
	for i=1:2
		if m.Amp[i] != 0.0 # Set stroke/hinge amplitude https://github.com/avikde/robobee3d/pull/127
			tcoord(i) .*= m.Amp[i]/currentAmpl(i) # scale pos
			tcoord(i+2) .*= m.Amp[i]/currentAmpl(i) # scale vel
		end
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
	if verbose
		println("Avg lift initial [mN]=", round(avgLift0, digits=4), ", in [mg]=", round(avgLift0*1000/9.81, digits=1))
	end
	return N, trajt, traj0, opt, avgLift0
end

"""Generate plot like in [Jafferis (2016)] Fig. 4"""
function openLoopPlot(m, opt, param0; save=false)
	function getResp(f, uamp, nlt=false)
		param = copy(param0)
		param[5] = nlt ? 2*param[2] : 0.0
		return createInitialTraj(m, opt, 0, f, [1e3, 1e2], param, 0; uampl=uamp, trajstats=true, thcoeff=0.1)
	end
	fs = 0.03:0.005:0.25
	mN_PER_V = 75/160

	p1 = plot(ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.8))
	p2 = plot(xlabel="Freq [kHz]", ylabel="Hinge ampl [deg]", legend=false, ylims=(0,100))

	function plotForTrans(nlt)
		nltstr = nlt ? "N" : "L"
		for Vamp=130:30:210
			println("Openloop @ ", Vamp, "V ", nltstr)
			uamp = Vamp*mN_PER_V
			amps = hcat(getResp.(fs, uamp, nlt)...)
			amps *= 180/pi # to degrees
			amps[1,:] /= (Vamp) # normalize
			amps[2,:] /= 2.0 # hinge ampl one direction
			# println(amps)
			plot!(p1, fs, amps[1,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
			plot!(p2, fs, amps[2,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
		end
	end

	plotForTrans(false)
	plotForTrans(true)

	plot(p1, p2, layout=(2,1), size=(500,500), dpi=200)
	if save
		savefig("olplot.png")
	end
	gui()
	error("Open loop plot")
end

"""Linear approx of wing AR constraint at cbar https://github.com/avikde/robobee3d/issues/113. Returns a,b s.t. dot(a, [cb,Aw]) <= b is the constraint"""
function wingARconstraintLin(cbar; maxAR=6)
	# AR<=?: see https://github.com/avikde/robobee3d/pull/105#issuecomment-562761586 originally
	# Note initial AR from Jafferis 2016 = 17/(54.59/17) ~= 5.3. should try to keep this
	fwingAR(cbAw::Vector) = -maxAR*cbAw[1]^2+cbAw[2]
	ptang(cb) = [cb, -fwingAR([cb, 0])]
	p0 = ptang(cbar)
	Dfwing(x) = ForwardDiff.gradient(fwingAR, x)
	return Dfwing(p0), dot(Dfwing(p0),p0)
end

"""Min lift involves Aw, cbar, dt. For now approx by fixing wing AR https://github.com/avikde/robobee3d/pull/119#issuecomment-577253382"""
function minLiftConstraintLin(minlift, param0, avgLift0)
	# Lift ~ (Aw/dt)^2 => calculate k needed
	cbar, τ1, mwing, wΨ, τ2, Aw, dt = param0
	Aw_dtmin = (Aw/dt)*sqrt(minlift/avgLift0)
	# constraint is already linear: -Aw + Aw_dtmin*dt <= 0
	return [-1.0, Aw_dtmin], 0
end

function trajMechPow(m, opt, traj, param)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
	mechpow = zeros(N)
	for k=1:N
		T = cu.transmission(m, traj[liy[:,k]], param; o2a=true)[2]
		mechpow[k] = traj[liy[3,k]]/T * traj[liu[1,k]]
	end
	return mechpow
end

"""One-off ID or opt"""
function opt1(m, traj, param, mode, minal, τ21ratiolim=2.0; testAffine=false, testAfter=false, testReconstruction=false, max_iter=4000, print_level=1, Φ=nothing)
	# A polytope constraint for the params: cbar >= cbarmin => -cbar <= -cbarmin. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	print(mode==2 ? "ID" : "Opt", " Φ=", isnothing(Φ) ? "-" : Φ, ", minal=", minal, ", τ2/1 lim=", τ21ratiolim, " => ")

    # cbar, τ1, mwing, wΨ, τ2, Aw, dt  = param
	# Poly constraint
	rholims = estimateWingDensity()
	# Could pick Wing AR constrain lin based on minal
	wARconstraintLinCbar = 2.0 + 0.5*minal # heuristic approx
	wARa, wARb = wingARconstraintLin(wARconstraintLinCbar; maxAR=4)
	mla, mlb = minLiftConstraintLin(minal, param0, avgLift0)
	# Polytope constraint
	Cp = Float64[0  0  0  0  0  mla[1]  mla[2]; # min lift
		0  -τ21ratiolim  0  0  1  0  0; # transmission nonlinearity τ2 <= τ21ratiolim * τ1
		0   0  -1  0  0  rholims[1]  0; # wing density mw >= Aw*ρ1
		0   0  1  0  0  -rholims[2]  0; # wing density mw <= Aw*ρ2
		wARa[1]   0  0  0  0  wARa[2]  0] # wing AR
	dp = [mlb; 0; 0; 0; wARb]

	if !isnothing(Φ)
		m.Amp[1] = deg2rad(Φ)
		# get new input traj
		traj = initTraj(m, param, KINTYPE; uampl=75, verbose=false)[3]
	end

	ret = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, Cp=Cp, dp=dp, print_level=print_level, max_iter=max_iter, testTrajReconstruction=testReconstruction)
	# append unactErr
	ret["unactErr"] = ret["eval_g"](ret["x"])[1:N] # 1 unact DOF
	ret["al"] = avgLift(m, opt, ret["traj"], ret["param"])
	uu = ret["traj"][(N+1)*ny:end]
	ret["u∞"] = norm(uu, Inf)
	ret["δact"] = maximum(abs.(actAng(m, opt, ret["traj"], ret["param"])))
	if testAfter
		cu.affineTest(m, opt, ret["traj"], ret["param"], POPTS)
	end
	# Calculate mechanical power
	ret["mechPow"] = trajMechPow(m, opt, ret["traj"], ret["param"])
	ret["comps"] = getComponents(m, opt, ret["traj"], ret["param"])
	ret["FD∞"] = norm(ret["comps"][6][1,:], Inf) # Also store drag (should be same as uinf for scaling but dynamics)
	
	println(ret["status"], ", ", round.(ret["param"]', digits=3), 
	", fHz=", round(1000/(N*ret["param"][end]), digits=1), 
	", al[mg]=", round(ret["al"] * 1000/9.81, digits=1), 
	", u∞=", round(ret["u∞"], digits=1), 
	", FD∞=", round(ret["FD∞"], digits=1), 
	", pow=", round(mean(abs.(ret["mechPow"])), digits=1), 
	", J=", round(ret["eval_f"](ret["x"]), digits=1), 
	", AR=", round(ret["param"][Aw_idx]/ret["param"][cb_idx]^2, digits=1), 
	", x=", round(ret["param"][Aw_idx]/ret["param"][cb_idx]*deg2rad(isnothing(Φ) ? 90 : Φ), digits=1))
	return ret
end

listOfParamTraj(rets...) = [ret["param"] for ret in rets], [ret["traj"] for ret in rets]

function getComponents(m::Wing2DOFModel, opt, traj1, param1)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj1)

	# Get the components
	yo, HMnc, HMc, HC, Hg, Hgact, HF, Hdamp, Hvel = cu.paramAffine(m, opt, traj1, param1, POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, param1)
	inertial = zeros(2,N)
	inertialc = similar(inertial)
	coriolis = similar(inertial)
	stiffdamp = similar(inertial)
	stiffdampa = similar(inertial)
	mechpow = zeros(N) # calculate from affine
	aero = similar(inertial)
	# to fill in for (non) inertial half of H https://github.com/avikde/robobee3d/pull/102
	npt1 = length(pt0)÷3
	H0 = zeros(size(inertial,1),npt1)

	for k=1:N
		y = yo(k)
		yn = yo(k+1)
		Ht(Hi) = cu.Hτ(Hi, y[1])
		# Divided up like this https://github.com/avikde/robobee3d/pull/119#issuecomment-577350049
		inertial[:,k] = [Ht(HMnc(y,yn) - HMnc(y,y))   H0  H0] * pt0
		inertialc[:,k] = [Ht(HMc(y,yn) - HMc(y,y) + HC(y))   H0  H0] * pt0
		stiffdamp[:,k] = [H0  Ht(Hdamp(y))  Ht(Hg(y))] * pt0
		stiffdampa[:,k] = [H0  H0  Ht(Hgact(y))] * pt0
		coriolis[:,k] = [Ht(HC(y))  H0  H0] * pt0
		aero[:,k] = [Ht(HF(y))  H0  H0] * pt0
		# put 
		Hh = [Ht(HMnc(y,yn) + HMc(y,yn) - HMnc(y,y) - HMc(y,y) + HC(y) + HF(y))   Ht(Hdamp(y))   Ht(Hg(y) + Hgact(y))]
		mechpow[k] = dot([H0[1,:]  Ht(Hvel(y))'  H0[1,:]], pt0) * dot(Hh[1,:], pt0)
	end
	return inertial, inertialc, coriolis, stiffdamp, stiffdampa, aero, mechpow 
end

"""Debug components in a traj"""
function debugComponentsPlot(m::Wing2DOFModel, opt, POPTS, ret)
	traj1, param1 = ret["traj"], ret["param"]
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, ret["traj"])
	inertial, inertialc, coriolis, stiffdamp, stiffdampa, aero, mechpow = ret["comps"]

	dt = param1[end]
	t2 = collect(0:(N-1))*dt

	function plotComponents(i, ylbl)
		# grav+inertial -- ideally they would "cancel" at "resonance"
		inertiastiff = inertial[i,:]+inertialc[i,:]+stiffdamp[i,:]+stiffdampa[i,:]
		tot = inertiastiff+aero[i,:]

		pl = plot(t2, (inertial[i,:] + inertialc[i,:]), linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, t2, stiffdamp[i,:], linewidth=2, label="g")
		# plot!(pl, t2, aero[i,:], linewidth=2, label="a")
		plot!(pl, t2, tot, linewidth=2, linestyle=:dash, label="act") # checked; this matches the row above ^
		if i==1
			plot!(pl, t2, stiffdampa[i,:], linewidth=2, label="ga")
			plot!(pl, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="")
		end

		pl2 = plot(t2, aero[i,:], linewidth=2, label="-dr", legend=:outertopright)
		plot!(pl2, t2, traj1[(N+1)*ny+1:end], linewidth=2, label="act")
		plot!(pl2, t2, coriolis[i,:], linewidth=2, label="cor")
		plot!(pl2, t2, inertiastiff, linewidth=2, label="is")
		hline!(pl2, [-1,1]*ret["FD∞"], ls=:dash, color=:black)
		
		pl3 = plot(t2, inertial[i,:], linewidth=2, label="inc", legend=:outertopright)
		plot!(pl3, t2, inertialc[i,:], linewidth=2, label="ic")
		plot!(pl3, t2, (inertial[i,:] + inertialc[i,:]), linewidth=2, linestyle=:dash, label="itot")
		plot!(pl3, t2, -(stiffdamp[i,:] + stiffdampa[i,:]), linewidth=2, label="-gtot")

		return pl, pl2, pl3
	end

	pl1 = plotTrajs(m, opt, listOfParamTraj(ret)...)
	pls, plcomp, plis = plotComponents(1, "stroke")
	plh, _, plih = plotComponents(2, "hinge")
	plot!(pl1[4], t2, 10*mechpow, label="mp", lw=2)
	plot!(pl1[4], t2, 10*ret["mechPow"], lw=2, ls=:dash, label="mpT")

	# Note that gamma is here
	# println("param = ", param1', ", Iw = ", param1[3] * (0.5 * param1[1])^2)
	return pl1[[1,2,4,6]]..., pls, plh, plcomp, plis, plih
end
