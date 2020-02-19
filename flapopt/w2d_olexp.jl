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
function statsFromAmplitudes(m, opt, Amp, param, Aw, Lw, freqHz)
	m.Amp .= Amp
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
	return [avgLift(m, opt, traj, param1); 0] # TODO: power
end

"Estimate lift, power for given stroke/hinge"
function calculateStats(m, opt, param, fname, Aw, Lw; spar0=[45., 45.])
	dat = readdlm(fname, ',', Float64, skipstart=1)
	spars = dat[:,4:5]
	# bitarray of rows where the spar proj angle has been recorded
	recordedpts = .!(spars[:,1] .≈ 0)
	# get mean wing pitch angle from the two spar measurements
	pitches = mean(hcat([wingPitchFromSparProjAng.(spars[recordedpts,i], Ref(spar0[i])) for i=1:2]...); dims=2)

	freqs = dat[recordedpts,2]
	strokes = dat[recordedpts,3]
	Amps = deg2rad.(hcat(strokes, 2*pitches))

	lifts = vcat([statsFromAmplitudes(m, opt, Amps[i,:], param, Aw, Lw, freqs[i]) for i =1:length(freqs)]'...)
	println(lifts)

	# # get new input traj
	# traj, Φ1 = initTraj(m, param, KINTYPE; uampl=75, verbose=false)[[3,5]]

	# println(statsFromAmplitudes(m, opt, Amp, param))
end

# --------------------------------------------------------

calculateStats(m, opt, param0, "data/normstroke/Param opt manuf 2 - halfbee1 a1.csv", 54, 15)

# # Main plot
# plot(
# 	olExpPlot2(
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")..., [120,160,190]), 
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [150,180,200]); 
# 		title="Wing 1A1"), 
# 	olExpPlot2(
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")..., [120,140,160]), 
# 		(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")..., [120,150,190]);
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

# # BigBee
# p = olExpPlot2(
# 	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee b1.csv")..., []), 
# 	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee 5b1.csv")..., []), 
# 	(readOLExpCSV("data/normstroke/Param opt manuf 2 - bigbee orig.csv")..., []); 
# 	title="BigBee", ulim=0.8)
# plot(p, size=(400,400))

# gui()
	