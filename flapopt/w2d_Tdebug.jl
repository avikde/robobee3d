using MAT

"""Use an open-loop sim to test for the nonlin transmission by producing normalized stroke, act disp, as well as lift/power dots.
WARNING: this overwrites POPTS"""
function openLoopTestTransmission(m, opt, param0)

	# # limit different params for component-wise optimization debugging
	# POPTS.plimsL .= copy(param0)
	# POPTS.plimsU .= copy(param0)
	# # Allow change in T1, Aw, mw
	# POPTS.plimsL[]

	"Returns strokeAmp, pitchAmp, actDisp, lift, power, "
	function getResp(f, uamp, nlt, lin, powplot=false)
		param = copy(param0)
		if nlt==1
			# this results in uinf = 73 but has nonlin transmission
			# ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 150; Φ=90, Qdt=3e4)
			param = [17.117, 2.469, 0.637, 4.477, 4.929, 70.894, 0.07] # 130, phi75
		elseif nlt==2
			param = [13.881, 2.469, 0.706, 2.867, 4.939, 78.581, 0.089] #1e3
		elseif nlt==3
			param = [12.637, 2.469, 0.632, 3.518, 4.939, 57.611, 0.068]
		elseif nlt==4
			param = [11.599, 2.469, 0.627, 3.348, 4.937, 54.311, 0.066]
			# param[3] = 0.75 # mw
			# param[2] *= 0.9
			# param[5] = 2*param[2]
		end
		if lin
			param[5] = 0
		end
		# println(param0, " ", param)
		ts = createInitialTraj(m, opt, 80, f, [1e3, 1e2], param, 212; uampl=uamp, trajstats2=true, h3=0.1, powplot=powplot)
		# To test, also use the amplitudes (same process as the lift/power estimates from data)
		Aw = param[6]
		cbar2 = param[1]
		Lw = Aw / sqrt(cbar2)
		return [ts; statsFromAmplitudes(m, opt, param, rad2deg.(ts[1:2]), param[6], Lw, f*1e3, uamp)]
	end
	fs = 0.03:0.01:0.22
	mN_PER_V = 75/180
	rightplot = false

	p1 = plot(ylabel=rightplot ? "" : "Norm. stroke ampl [deg/V]", ylims=(0.3,0.8), legend=:topleft)
	p2 = plot(xlabel="Freq [kHz]", ylabel="Pitch ampl [deg]", legend=false, ylims=(0,80))
	p3 = plot(ylabel=rightplot ? "" : "Act. disp [mm]", legend=false, xlabel="Freq [kHz]")
	# if rightplot
	# 	yaxis!(p1, false)
	# 	yaxis!(p3, false)
	# end
	p4 = plot(xlabel="Lift [mg]", ylabel="Power [mW]", legend=false)
	p5 = plot(xlabel="lift est", ylabel="power est")
	p6 = plot(xlabel="Freq [kHz]", ylabel="Lift [mg]", legend=false)
	p7 = plot(xlabel="Freq [kHz]", ylabel="Power [mW]", legend=false)
	mss = [:circle, :utriangle, :dtriangle, :star, :rect]

	function plotForTrans(nlt, Vamp, lin=false; powplot=false)
		nltstr = nlt == 0 ? "Orig" : string(nlt, lin ? "L" : "N")
		actdisps = Dict{Float64, Float64}()
		ms = mss[nlt+1]
		lbl = string(nltstr, " ", Vamp,"V")
		
		println("Openloop @ ", Vamp, "V ", nltstr)
		uamp = Vamp*mN_PER_V
		amps = hcat(getResp.(fs, uamp, nlt, lin, powplot)...)
		amps[1:2,:] *= 180/pi # to degrees
		amps[1,:] /= (Vamp) # normalize
		amps[2,:] /= 2.0 # hinge ampl one direction
		# println(amps)
		plot!(p1, fs, amps[1,:], lw=2, label=lbl, markershape=ms)
		plot!(p2, fs, amps[2,:], lw=2, label=lbl, markershape=ms)
		plot!(p3, fs, amps[3,:], lw=2, label=lbl, markershape=ms)
		plot!(p6, fs, amps[4,:], lw=2, label=lbl, markershape=ms)
		plot!(p7, fs, amps[5,:], lw=2, label=lbl, markershape=ms)
		# pick the one that produced the max lift
		imax = argmax(amps[4,:])
		scatter!(p4, [amps[4,imax]], [amps[5,imax]], label=lbl, markershape=ms, markersize=(nlt==0 || lin) ? 4 : 8)
		# println(nlt, " max at f=", fs[imax], ", Vamp=", Vamp, " ", imax, " ", amps[4,:])
		# scatter!(p5, [amps[6,imax]], [amps[7,imax]], label=lbl, markershape=ms)
	end

	plotForTrans(0, 178)
	plotForTrans(1, 230)
	# plotForTrans(2, 215)
	plotForTrans(3, 190)
	plotForTrans(4, 183)
	plotForTrans(1, 215, true)
	# plotForTrans(2, 200, true)
	plotForTrans(3, 180, true)
	plotForTrans(4, 172, true)

	# return plot(p1, p2, p3, p4, size=(600,500))
	return plot(p1, p2, p3, p4, p6,p7, size=(600,500))
end
openLoopTestTransmission(mop...)

## New idea: --------------------------------------------------
#  (a) select some weights and
# run opt. plot some components (maybe power, speed, Fact),
# (b) run the same opt with a tau2 ratio lim to force linear,
# (c) same params as (a) but do a sweep with the same peak
# voltage to obtain the best results with linear transmission

# include("run_w2d.jl") # just once
# includet("w2d_paramopt.jl")

function Tdebug(ret1, minal, Φ, Qdt; fullDicts=false)
	"Sweep over freq and pick the one with the best lift"
	function sweepFreqAndPickBestLift(param, uamp)
		ts2arr = hcat([[f; createInitialTraj(m, opt, 80, f, [1e3, 1e2], param, 212; uampl=uamp, trajstats2=true, h3=0.1)] for f in range(0.15,0.2,step=0.1)]...)
		# ^ rows of returned mat are [freq; stroke ampl; pitch ampl; act disp; avglift; pow]
		# change row 5 into sp. lift
		ts2arr[5,:] ./= (ts2arr[4,:] * uamp)
		imax = argmax(ts2arr[5,:])
		println("Linear best was at ", round.(ts2arr[:,imax]', digits=2))
		return ts2arr[5:6,imax] # sp. lift, pow
	end
	# (a) nonlin
	retNL = @time opt1(m, ret1["traj"], ret1["param"], 1, minal, 2.0; Φ=Φ, Qdt=Qdt)
	# (b) lin opt
	retL = @time opt1(m, ret1["traj"], ret1["param"], 1, minal, 0; Φ=Φ, Qdt=Qdt)
	# # (c) lin at same params as (a)
	# paramL2 = copy(retNL["param"])
	# L2 = sweepFreqAndPickBestLift(paramL2, retNL["u∞"])
	# println(retNL["param"]', retL["param"]')

	if fullDicts
		return retNL, retL
	end

	mp(ret) = mean(cu.ramp.(ret["mechPow"]))
	mact(ret) = ret["u∞"] * ret["δact"]

	return vcat([Qdt, retNL["al"]/mact(retNL),  retL["al"]/mact(retL), mp(retNL),  mp(retL)])#, L2)
end

function TdebugTraj(pstroke, pact, pforce, ppow, plift, pdrag, ret, lbl)
	println("Params = ", round.(ret["param"]',digits=3))
	inertial, inertialc, coriolis, stiffdamp, stiffdampa, aero, mechpow = ret["comps"]
	actt = actAng(m, opt, ret["traj"], ret["param"])
	liftt = trajAero(m, opt, ret["traj"], ret["param"], :lift)
	dt = ret["param"][end]
	t = 0:dt:(N)*dt
	t2 = collect(0:(N-1))*dt
	plot!(pact, t, actt, label=lbl, lw=2)
	plot!(pstroke, t, ret["traj"][1:ny:(N+1)*ny], label=lbl, lw=2)
	plot!(pforce, t2, ret["traj"][(N+1)*ny+1:end], label=lbl, lw=2)
	plot!(ppow, t2, ret["mechPow"], label=lbl, lw=2)
	plot!(plift, t2, liftt, label=lbl, lw=2)
	plot!(pdrag, t2, trajAero(m, opt, ret["traj"], ret["param"], :drag), label=lbl, lw=2)
end
##

res = hcat([Tdebug(ret1, 130, nothing, Qdt) for Qdt in range(1e4,1e5,length=8)]...)
matwrite("Tdebug1.zip", Dict("results"=>res); compress=true)
##
retNL, retL = Tdebug(ret1, 130, nothing, 1e4; fullDicts=true)
retNL2, retL2 = Tdebug(ret1, 130, nothing, 1e5; fullDicts=true)
## 

p1 = plot(res[1,:], res[2,:], xaxis=:log,  xlabel="Weighting on power", ylabel="Lift/mact [ ]", label="NL", lw=2, markershape=:circle, legend=false)
plot!(p1, res[1,:], res[3,:], label="L", lw=2, markershape=:utriangle)
# plot!(p1, res[1,:], res[6,:], label="L2", lw=2, markershape=:dtriangle)

p2 = plot(res[1,:], res[4,:], xaxis=:log, xlabel="Weighting on power", ylabel="Avg power [mW]", label="NL", lw=2, markershape=:circle, legend=:topleft)
plot!(p2, res[1,:], res[5,:], label="L", lw=2, markershape=:utriangle)
# plot!(p2, res[1,:], res[7,:], label="L2", lw=2, markershape=:dtriangle)
#
pstroke = plot(ylabel="Stroke [rad]", legend=:bottomright)
pact = plot(ylabel="Act disp [mm]", legend=false)
pforce = plot(ylabel="Act force [mN]", legend=false)
ppow = plot(ylabel="Pow [mW]", legend=false)
plift = plot(xlabel="t [ms]", ylabel="Lift [mg]", legend=false)
pdrag = plot(xlabel="t [ms]", ylabel="Drag [mg]", legend=false)
TdebugTraj(pstroke, pact, pforce, ppow, plift, pdrag, retNL, "NL1")
TdebugTraj(pstroke, pact, pforce, ppow, plift, pdrag, retL, "L1")
TdebugTraj(pstroke, pact, pforce, ppow, plift, pdrag, retNL2, "NL2")
TdebugTraj(pstroke, pact, pforce, ppow, plift, pdrag, retL2, "L2")

plot(p1, p2, pstroke, pact, pforce, ppow, plift, pdrag, layout=(4,2), size=(500,800))
# plot(p1, pact, p2, pforce)
# pl1 = debugComponentsPlot(m, opt, POPTS, retNL)
# plot(pl1..., size=(800,600))
# pl = scatter(xlabel="Lift/mact", ylabel="Pow/mact")
# TdebugPlot!(pl, Tdebug(ret1, 130, 75, 3e4))
# pls = debugComponentsPlot(m, opt, POPTS, ret2)
# plot(pls..., size=(800,600))
# gui()
