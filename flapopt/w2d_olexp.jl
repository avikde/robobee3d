using Plots, DelimitedFiles, Statistics

include("w2d_model.jl")
include("w2d_paramopt.jl") # for calculating stats for measurements

# Simulation traces for norm stroke ------------------------------------

"Normalize by lowering T1 corresponding to the maximum output disp"
function nlNormalizeByOutput(m, opt, param, σa0, T1scale; τ2ratio=2.0)
	if !isnothing(T1scale)
		param[2] *= T1scale
	else
		# Automatically try to match output. From transmission:
		# τfun = σa -> τ1*σa + τ2/3*σa^3
		# = τ1*σa + τ1*τ2ratio/3*σa^3 = τ1*σa*(1 + τ2ratio/3*σa^2)
		# So scale τ1 by 1/(1 + τ2ratio/3*σa^2)
		param[2] *= 1/(1 + τ2ratio/3*σa0^2)
	end
	param[5] = τ2ratio*param[2]
	return param
end

function simNormStroke!(pls, plh, pla, mop, param, frange, Vs; mN_PER_V=75/160, kwargs...)
	function simVf(fHz, Vamp)
		amps = createInitialTraj(mop[1:2]..., 0, fHz*1e-3, [1e3, 1e2], param, 0; uampl=75, trajstats=true, h3=0.1)
		amps[1:2] *= 180/pi # to degrees
		amps[1] /= (Vamp) # normalize
		amps[3] *= (1000/Vamp)
		amps[2] /= 2.0 # hinge ampl one direction
		return amps
	end

	# frequency sweep
	fs = range(frange..., length=20)
	for Vamp in Vs
		println("Running V=", round(Vamp))
		amps = hcat(simVf.(fs, Ref(Vamp))...)
		
		plot!(pls, fs, amps[1,:], lw=2, label=string(Vamp,"V"); kwargs...)
		if !isnothing(plh)
			plot!(plh, fs, amps[2,:], lw=2, label=string(Vamp,"V"); kwargs...)
		end
		if !isnothing(pla)
			plot!(pla, fs, amps[3,:], lw=2, label=string(Vamp,"V"); kwargs...)
		end
	end
end

"""Generate plot like in [Jafferis (2016)] Fig. 4"""
function openLoopPlot(m, opt, param0, Vamps; save=false, NLT1scales=nothing, rightplot=false)
	function getResp(f, uamp, nlt, σa0=nothing, T1scale=nothing)
		param = copy(param0)
		if nlt==1
			param = nlNormalizeByOutput(m, opt, param, σa0, T1scale)
		elseif nlt==2
			param[POPTS.τinds[1]] *= T1scale # and don't change \tau2
		end
		ts = createInitialTraj(m, opt, 0, f, [1e3, 1e2], param, 0; uampl=uamp, trajstats=true, h3=0.1)
		# println("act disp=",ts[end])
		return ts
	end
	fs = 0.03:0.005:0.25
	mN_PER_V = 75/160

	p1 = plot(ylabel=rightplot ? "" : "Norm. stroke ampl [deg/V]", ylims=(0.3,0.8), legend=false, title=rightplot ? "High inertia" : "Low inertia")
	p2 = plot(xlabel="Freq [kHz]", ylabel="Hinge ampl [deg]", legend=false, ylims=(0,100))
	p3 = plot(ylabel=rightplot ? "" : "Norm. act. disp [um/V]", legend=rightplot ? false : :topleft, ylim=(1.0,2.4), xlabel="Freq [kHz]")
	# if rightplot
	# 	yaxis!(p1, false)
	# 	yaxis!(p3, false)
	# end
	function plotForTrans(nlt; T1scales=nothing)
		nltstr = nlt == 1 ? "N" : (nlt == 2 ? "LL" : "L")
		actdisps = Dict{Float64, Float64}()
		
		for ii=1:length(Vamps)
			Vamp = Vamps[ii]
			T1scale = isnothing(T1scales) ? nothing : T1scales[ii]

			print("Openloop @ ", Vamp, "V ", nltstr)
			uamp = Vamp*mN_PER_V
			σa0 = nothing
			# println("HI", σa0)
			amps = hcat(getResp.(fs, uamp, nlt, σa0, T1scale)...)
			
			actdisps[Vamp] = maximum(amps[3,:])
			println(", act disp=", round(actdisps[Vamp], digits=3), "mm")
			amps[1:2,:] *= 180/pi # to degrees
			amps[1,:] /= (Vamp) # normalize
			amps[3,:] *= (1000/Vamp)
			amps[2,:] /= 2.0 # hinge ampl one direction
			# println(amps)
			plot!(p1, fs, amps[1,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt==1 ? :solid : (nlt == 2 ? :dot : :dash))
			plot!(p2, fs, amps[2,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt==1 ? :solid : (nlt == 2 ? :dot : :dash))
			plot!(p3, fs, amps[3,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt==1 ? :solid : (nlt == 2 ? :dot : :dash))
		end
		return actdisps
	end

	plotForTrans(0)
	plotForTrans(1; T1scales=NLT1scales)
	# plotForTrans(2; T1scales=NLT1scales) # NOTE: for Rob test

	println("dens=", param0[3]/param0[6], ", koratio=", m.kbo[1]/(m.kbo[1] + m.ka/param0[2]^2))
	# plot(p1, #= p2, layout=(2,1),  =# size=(400, 300), dpi=200)
	return plot(p1, p3, layout=(2,1), size=(400,250))
end

function openLoopPlotFinal(m, opt, param0)
	# mw=0.55 -> [0.975, 0.925], mw=0.7 -> [0.96, 0.9]
	param = copy(param0)
	param[3] = 0.55 #mw
	pl1 = openLoopPlot(m, opt, param, range(140, 180,length=2); NLT1scales=[0.975, 0.925])
	param[3] = 0.7 #mw
	pl2 = openLoopPlot(m, opt, param, range(140, 180,length=2); NLT1scales=[0.96, 0.9], rightplot=true)
	plot(pl1, pl2, size=(500, 500), dpi=150)
	gui()
	error("Open loop plot")
end

# Experimental data ------------------------------------

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

STROKE_ACT_COORDS = false # to respond to Rob question about what act disp was in these trials

function toAct(stroke, τcoeffs)
	param = copy(param0)
	param[POPTS.τinds[1]] = τcoeffs[1]
	param[POPTS.τinds[2]] = τcoeffs[2]
	return [cu.transmission(m, [deg2rad(stroke[k]/2), 0, 0, 0], param; o2a=true)[1][1] for k=1:length(stroke)]
end

function olExpPlotCurves!(p, dataset, lbl; kwargs...)
	V, stroke, wingDims, showV, ms, τcoeffs = dataset
	Nv = length(V)
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], (STROKE_ACT_COORDS ? toAct(stroke[i][:,2], τcoeffs) : stroke[i][:,2]/V[i]), label=string(lbl,Int(V[i])), markershape=ms, lw=2.5; kwargs...)
			# # Plot sim data overlaid
			# plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=string(lbl,Int(V[i])), markershape=ms, lt=:scatter; kwargs...)
			# param = copy(mop[end])
			# param[6] = wingDims[1]
			# param[1] = (wingDims[1]/wingDims[2])^2
			# simNormStroke!(p, nothing, nothing, mop, param, (120,200), [V[i]]; legend=false, kwargs...)
		end
	end
end

function olExpPlot2(mop, datasets...; title="", ulim=0.5)
	p = STROKE_ACT_COORDS ? 
		plot(xlabel="Freq [Hz]", ylabel="Peak act disp [mm]", legend=:topleft, title=title) :
		plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,ulim), legend=:topleft, title=title)

	@assert length(datasets) > 0
	lss = [:solid, :dash, :dot, :dash]
	for k=1:length(datasets)
		olExpPlotCurves!(p, datasets[k], "", ls=lss[k])
	end

	# # TEST MODEL FIT
	# param = [3.2^2,  # cbar2[mm^2] (area/R)^2
	# 		2, # τ1 (from 3333 rad/m, [Jafferis (2016)])
	# 		0.8, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
	# 		1.5, # wΨ [mm]
	# 		0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
	# 		70, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
	# 		0.0758 # dt DOES NOT AFFECT
	# 	]
	# simNormStroke!(p, nothing, nothing, mop, param, (100,200), [160]; mN_PER_V=75/180, legend=false)

	return p
end

"Get the pitch angle from this measurement: see https://github.com/avikde/robobee3d/issues/145#issuecomment-588361995"
wingPitchFromSparProjAng(angProj, ang0=45) = rad2deg(acos(1 - tan(deg2rad(angProj)) / tan(deg2rad(ang0))))

"Only Aw, cbar (wing dims), and dt (freq) are use to calculate lift. Use param0 otherwise."
function statsFromAmplitudes(m, opt, param, tau, Amp, Aw, Lw, freqHz, Fact)
	m.Amp .= Amp
	# m.Amp[1] = Amp[1]
	N, traj = initTraj(m, param, 1; uampl=75, verbose=false)[[1,3]]
	# scale velocities according to the new freq
	forig = 1000/(N*param[end])
	ny = 4
	traj[3:ny:(N+1)*ny] .*= freqHz/forig
	traj[4:ny:(N+1)*ny] .*= freqHz/forig
	# update the params according to the inputs
	param1 = copy(param)
	param1[1] = (Aw/Lw)^2
	param1[6] = Aw
	param1[end] = 1000/(N*freqHz)

	# pls = plotTrajs(m, opt, [param], [traj])
	# plot(pls...)
	# gui()

	# FD0 = maximum(filter(!isnan, trajAero(m, opt, traj, param1, :drag)))
	# Approximate power with least dependence on dynamical parameters. If this and force were in-phase and sinusoidal, then
	# Pmech ~ Fact/T*Lw*Phi*f/pi (1/pi = integral of sin^2 x)
	Fwing = Fact/tau
	println(Fact, "/", tau, " * ", Amp[1], " * 1e-3*", freqHz, " = ", Fwing*Amp[1]*(freqHz*1e-3))
	return [avgLift(m, opt, traj, param1); Fwing*Amp[1]*(freqHz*1e-3)]
end

"Estimate lift, power for given stroke/hinge"
function calculateStats(tau, mop, fname, Aw, Lw, spar10, spar20; kVF=75/180)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	spars = dat[:,4:5]
	# bitarray of rows where the spar proj angle has been recorded
	recordedpts = .!(spars[:,1] .≈ 0)
	spar0 = (spar10, spar20)
	# get mean wing pitch angle from the two spar measurements
	pitches = mean(hcat([wingPitchFromSparProjAng.(spars[recordedpts,i], Ref(spar0[i])) for i=1:2]...); dims=2)

	volts = dat[recordedpts,1]
	freqs = dat[recordedpts,2]
	strokes = dat[recordedpts,3]
	# Here Amps = [stroke amplitude rad,  2 * pitch amplitude rad] (stacked column)
	Amps = deg2rad.(hcat(strokes, 2*pitches))

	println("file: ", fname)
	stats = vcat([
		[volts[i]; freqs[i]; statsFromAmplitudes(mop..., tau, Amps[i,:], Aw, Lw, freqs[i], volts[i])]'
		for i = 1:length(freqs)]...)
	
	println(fname, round.(pitches, digits=1), round.(stats[:,3], digits=1), round.(stats[:,4], digits=1))
	return stats, volts
end

"Tuples of Aw, Lw, spar1, spar2 from https://github.com/avikde/robobee3d/issues/145. Measured by roughly using area tool in Autocad, and measuring the longest membrane length"
bbHingeOffs = 1.3
sdabHingeOffs = 2.3
wingRootAddedOffset = 3.28 # to get to the R=17 in Jafferis et al 2016
wingDims = Dict{String,Tuple{Float64,Float64,Float64,Float64}}(
	"1a" => (54.4, 12.42 + sdabHingeOffs + wingRootAddedOffset, 45.45, 45.8),
	"1al" => (121.76, 18.8 + sdabHingeOffs + wingRootAddedOffset, 45.45, 45.8),
	# "1alx2" => (243.52, 18.8 + sdabHingeOffs + wingRootAddedOffset, 45.45, 45.8),
	"1b" => (64.4, 14.51 + sdabHingeOffs + wingRootAddedOffset, 36.18, 45.8),
	"4b" => (54.4, 13.8 + sdabHingeOffs + wingRootAddedOffset, 45, 46),
	"4b2" => (54.4, 13.8 + sdabHingeOffs + wingRootAddedOffset, 45, 46),
	"5b" => (58, 15 + sdabHingeOffs + wingRootAddedOffset, 47, 48),
	"bigbee" => (179, 25.5 + bbHingeOffs + wingRootAddedOffset, 45, 45), # measured https://github.com/avikde/robobee3d/issues/145
	"4l" => (107.3, 18.65 + sdabHingeOffs + wingRootAddedOffset, 45, 46),
)
τCOEFFS_MOD1 = (1.8, 4.4)
τCOEFFS_SDAB = (2.6666, 0)
τCOEFFS_BB = (3.28, 0)

function liftPowerPlot(mop; includeBigbee=true)
	p1 = plot(xlabel="FL [mg]", ylabel="pow (S*V*f)", legend=:topleft, title="SDAB actuator", ylim=(5,25))
	p3 = plot(xlabel="FL/V^2 [ug/V^2]", ylabel="pow (S*V*f)", legend=false, ylim=(5,25))

	function addToPlot!(p, pn, lbl, tau, args...; kwargs...)
		s, Vs = calculateStats(tau, args...)
		# # Old: plot all voltages together
		# scatter!(p, s[:,3], s[:,4], label=lbl; ms=6, kwargs...)
		# scatter!(pn, 1000*s[:,3]./s[:,1].^2, s[:,4], label=lbl; ms=6, kwargs...)
		# New: group by voltage
		uVs = unique(Vs)
		for uV in uVs
			lbl2 = string(lbl, " ", convert(Int,uV),"V")
			inds = Vs .== uV
			scatter!(p, s[inds,3], s[inds,4], label=lbl2; ms=6, kwargs...)
			scatter!(pn, 1000*s[inds,3]./s[inds,1].^2, s[inds,4], label=lbl2; ms=6, kwargs...)
		end
	end
	# addToPlot!(p1, p3, "hb 1a1", mop, "data/normstroke/Param opt manuf 2 - halfbee1 a1.csv", wingDims["1a"]...)
	addToPlot!(p1, p3, "orig", τCOEFFS_SDAB[1], mop, "data/normstroke/Param opt manuf 2 - sdab1.csv", wingDims["1a"]...; markershape=:circle)

	τmod1Avg = 0.5 * (τCOEFFS_MOD1[1] + (τCOEFFS_MOD1[1] + τCOEFFS_MOD1[2] * 0.35^2))
	# addToPlot!(p1, p3, "sdab1", mop, "data/normstroke/Param opt manuf 2 - beckysdab.csv", wingDims["1a"]...; markershape=:circle)
	# addToPlot!(p1, p3, "modreg", mop, "data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv", wingDims["1a"]...; markershape=:rect)
	# ^ was originally included
	# addToPlot!(p1, p3, "mod1 4b1", mop,  "data/normstroke/Param opt manuf 2 - mod4 b h1.csv", wingDims["4b"]...; markershape=:dtriangle)
	addToPlot!(p1, p3, "modhAR", τmod1Avg, mop,  "data/normstroke/Param opt manuf 2 - mod4 b h2.csv", wingDims["4b"]...; markershape=:utriangle)

	if includeBigbee
		p2 = plot(xlabel="FL [mg]", ylabel="pow (S*V*f)", legend=:topleft, title="BigBee actuator", ylim=(2,15))
		p4 = plot(xlabel="FL/V^2 [ug/V^2]", ylabel="pow (S*V*f)", legend=false, ylim=(2,15))
		addToPlot!(p2, p4, "bigbee orig", τCOEFFS_BB[1], mop, "data/normstroke/Param opt manuf 2 - bigbee orig.csv", wingDims["bigbee"]...)
		addToPlot!(p2, p4, "bigbee 1b", τCOEFFS_BB[1], mop,  "data/normstroke/Param opt manuf 2 - bigbee b1.csv", wingDims["1b"]...; markershape=:rect)
		addToPlot!(p2, p4, "bigbee 4l", τCOEFFS_BB[1], mop, "data/normstroke/Param opt manuf 2 - bigbee 4l3.csv", wingDims["4l"]...; markershape=:utriangle)
		# addToPlot!(p2, p4, "bbx 1al", mop, "data/normstroke/Param opt manuf 2 - bbx 1alx2.csv", wingDims["1alx2"]...; markershape=:star5)
		# addToPlot!(p2, p4, "bbx o", mop,  "data/normstroke/Param opt manuf 2 - bbx bborig.csv", wingDims["bigbee"]...; markershape=:dtriangle)
		# addToPlot!(p2, p4, "bbx 4l", mop, "data/normstroke/Param opt manuf 2 - bbx 4l3.csv", wingDims["4l"]...; markershape=:dtriangle)
		addToPlot!(p2, p4, "bigbee oA", τCOEFFS_BB[1], mop,  "data/normstroke/Param opt manuf 2 - bigbee originalA.csv", wingDims["bigbee"]...; markershape=:dtriangle)
		# addToPlot!(p2, p4, "bigbee 1al", mop,  "data/normstroke/Param opt manuf 2 - bigbee 1al.csv", wingDims["1al"]...; markershape=:dtriangle)
		plot(p1, p3, p2, p4, size=(800,500))
	else
		plot(p1, p3, size=(800,250))
	end
end

normStrokeSDAB(mop) = plot(
	olExpPlot2(mop, 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")..., wingDims["1a"], [120,160,190], :rect, τCOEFFS_MOD1), 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - sdab1.csv")..., wingDims["1a"], [], :circle, τCOEFFS_SDAB); 
		# (readOLExpCSV("data/normstroke/Param opt manuf 2 - beckysdab.csv")..., [], :circle); 
		title="Wing 1A1"), 
	olExpPlot2(mop, 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h2.csv")..., wingDims["4b2"], [120,150,190], :utriangle, τCOEFFS_MOD1),
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")..., wingDims["4b"], [120,150,190], :+, τCOEFFS_SDAB),
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")..., wingDims["4b"], [120,140], :dtriangle, τCOEFFS_MOD1);
		title="Wing 4B1"),
	size=(800,400))

normStrokeBigBee(mop) = plot(
	olExpPlot2(mop, 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee b1.csv")..., wingDims["1b"], [], :utriangle, τCOEFFS_BB), 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee orig.csv")..., wingDims["bigbee"], [], :circle, τCOEFFS_BB), 
		# (readOLExpCSV("data/normstroke/Param opt manuf 2 - bbx bborig.csv")..., wingDims["bigbee"], [], :dtriangle, τCOEFFS_BB),  
		# (readOLExpCSV("data/normstroke/Param opt manuf 2 - bbx 1alx2.csv")..., wingDims["1alx2"], [], :dtriangle, τCOEFFS_BB), 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee 4l3.csv")..., wingDims["4l"], [150,180,200], :rect, τCOEFFS_BB), 
		(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee originalA.csv")..., wingDims["bigbee"], [], :star5, τCOEFFS_BB); 
		title="BigBee", ulim=0.75), 
	size=(600,400))
##

# --------------------------------------------------------
mop = (m, opt, param0)

liftPowerPlot(mop)

# normStrokeSDAB(mop)
# normStrokeBigBee(mop)

# openLoopPlotFinal(mop...)

gui()
