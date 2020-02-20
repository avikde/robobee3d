using Plots, DelimitedFiles, Statistics

include("w2d_model.jl")
include("w2d_paramopt.jl") # for calculating stats for measurements

# for comparison; my SDAB vs. others. these are all at 180
becky3L = [140	50.6;
	145	50.74;
	150	51.38;
	155	51.79;
	160	57.75;
	165	67.38;
	170	75.51;
	175	74.89;
	180	72.94]
becky3R = [140	52.88;
	145	54.21;
	150	56.47;
	155	59.6;
	160	69;
	165	73.6;
	170	72.65;
	175	72.46;
	180	68.01]
patrickL = [120	50.3;
	125 50.5;
	130 50;
	135 50.3;
	140	50.2;
	145 50.2;
	150	50.3;
	155	51.5;
	160	53.2;
	165 59;
	168 62;
	170	63.4]
patrickR = [120	47.7;
	125 47.4;
	130 46.8;
	135 46.8;
	140	49.1;
	145 49.4;
	150	49;
	155	53.7;
	160	59;
	165 65;
	168 67;
	170	68]

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlotCurves!(p, dataset, lbl; kwargs...)
	V, stroke, showV = dataset
	Nv = length(V)
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=string(lbl,Int(V[i])), markershape=:auto, lw=2.5; kwargs...)
		end
	end
end

function olExpPlot2(datasets...; title="", ulim=0.6)
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,ulim), legend=:topleft, title=title)

	@assert length(datasets) > 0
	olExpPlotCurves!(p, datasets[1], "")
	if length(datasets) > 1
		olExpPlotCurves!(p, datasets[2], ""; ls=:dash)
	end
	if length(datasets) > 2
		olExpPlotCurves!(p, datasets[3], ""; ls=:dot)
	end
	return p
end

"Get the pitch angle from this measurement: see https://github.com/avikde/robobee3d/issues/145#issuecomment-588361995"
wingPitchFromSparProjAng(angProj, ang0=45) = rad2deg(acos(1 - tan(deg2rad(angProj)) / tan(deg2rad(ang0))))

"Only Aw, cbar (wing dims), and dt (freq) are use to calculate lift. Use param0 otherwise."
function statsFromAmplitudes(m, opt, param, Amp, Aw, Lw, freqHz, Fact)
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
	# Pmech ~ Fact*Lw*Phi*f/pi (1/pi = integral of sin^2 x)
	return [avgLift(m, opt, traj, param1); Fact*Amp[1]*(freqHz*1e-3)]
end

"Estimate lift, power for given stroke/hinge"
function calculateStats(mop, fname, Aw, Lw, spar10, spar20; kVF=75/180)
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
	Amps = deg2rad.(hcat(strokes, 2*pitches))

	stats = vcat([
		[volts[i]; freqs[i]; statsFromAmplitudes(mop..., Amps[i,:], Aw, Lw, freqs[i], volts[i])]'
		for i = 1:length(freqs)]...)
	
	println(fname, round.(pitches, digits=1), round.(stats[:,3], digits=1), round.(stats[:,4], digits=1))
	return stats
end

"Tuples of Aw, Lw, spar1, spar2 from https://github.com/avikde/robobee3d/issues/145. Measured by roughly using area tool in Autocad, and measuring the longest membrane length"
wingDims = Dict{String,Tuple{Float64,Float64,Float64,Float64}}(
	"1a" => (54, 12.42, 45.45, 45.8),
	"1b" => (64.4, 14.51, 36.18, 45.8),
	"4b" => (54.2, 13.8, 45, 46),
	"5b" => (58, 15, 47, 48),
	"bigbee" => (156, 25.5, 45, 45),
	"4l" => (107.3, 18.65, 45, 46),
)

function liftPowerPlot(mop)
	p1 = plot(xlabel="FL [mg]", ylabel="pow (S*V*f)", legend=:bottomright, title="SDAB actuator")
	p2 = plot(xlabel="FL [mg]", ylabel="pow (S*V*f)", legend=:topright, title="BigBee actuator")
	p3 = plot(xlabel="FL/V^2 [ug/V^2]", ylabel="pow (S*V*f)", legend=false)
	p4 = plot(xlabel="FL/V^2 [ug/V^2]", ylabel="pow (S*V*f)", legend=false)

	function addToPlot!(p, pn, lbl, args...; kwargs...)
		s = calculateStats(args...)
		scatter!(p, s[:,3], s[:,4], label=lbl; ms=6, kwargs...)
		scatter!(pn, 1000*s[:,3]./s[:,1].^2, s[:,4], label=lbl; ms=6, kwargs...)
	end
	addToPlot!(p1, p3, "hb 1a1", mop, "data/normstroke/Param opt manuf 2 - halfbee1 a1.csv", wingDims["1a"]...)
	addToPlot!(p1, p3, "mod1 1a1", mop, "data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv", wingDims["1a"]...; markershape=:rect)
	addToPlot!(p1, p3, "mod1 4b1", mop,  "data/normstroke/Param opt manuf 2 - mod4 b h1.csv", wingDims["4b"]...; markershape=:utriangle)
	addToPlot!(p1, p3, "mod1 4b2", mop,  "data/normstroke/Param opt manuf 2 - mod4 b h2.csv", wingDims["4b"]...; markershape=:dtriangle)
	addToPlot!(p2, p4, "bigbee orig", mop, "data/normstroke/Param opt manuf 2 - bigbee orig.csv", wingDims["bigbee"]...)
	addToPlot!(p2, p4, "bigbee 4l", mop, "data/normstroke/Param opt manuf 2 - bigbee 4l3.csv", wingDims["4l"]...; markershape=:rect)
	addToPlot!(p2, p4, "bigbee 1b", mop,  "data/normstroke/Param opt manuf 2 - bigbee b1.csv", wingDims["1b"]...; markershape=:utriangle)

	plot(p1, p2, p3, p4, size=(600,600))
end

# --------------------------------------------------------
mop = (m, opt, param0)

# liftPowerPlot(mop)

# # Main plot
# plot(
# 	olExpPlot2(
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")..., [120,160,190]), 
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [150,180,200]); 
# 		title="Wing 1A1"), 
# 	olExpPlot2(
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h2.csv")..., [120,150,200]),
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")..., [120,150,190]),
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")..., [120,140,160]);
# 		title="Wing 4B1"),
# 	size=(800,400))

# # plot of comparing different SDAB
# p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.6), legend=:topleft, title="Different SDAB")
# olExpPlotCurves!(p, readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [180], "Avik ")
# olExpPlotCurves!(p, ([180.0], [becky3L], []), "Becky 3L ")
# olExpPlotCurves!(p, ([180.0], [becky3R], []), "Becky 3R ")
# olExpPlotCurves!(p, ([180.0], [patrickL], []), "Patrick L ")
# olExpPlotCurves!(p, ([180.0], [patrickR], []), "Patrick R ")
# plot(p, size=(400,400))

# BigBee
p = olExpPlot2(
	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee b1.csv")..., []), 
	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee 4l3.csv")..., [150,180,200]), 
	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee orig.csv")..., []); 
	title="BigBee", ulim=0.75)
plot(p, size=(400,400))

gui()
	