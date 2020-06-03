
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


