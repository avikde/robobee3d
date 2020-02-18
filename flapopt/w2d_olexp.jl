using Plots, DelimitedFiles

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
patrick = [120	47.7;
	140	49.1;
	150	49;
	155	53.8;
	160	59;
	170	68]

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlotCurves!(p, V, stroke, showV, lbl; kwargs...)
	Nv = length(V)
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=string(lbl,Int(V[i])), markershape=:auto, lw=2.5; kwargs...)
		end
	end
end

function olExpPlot(V1, stroke1, showV1, V2, stroke2, showV2; title="")
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.6), legend=:topleft, title=title)

	olExpPlotCurves!(p, V1, stroke1, showV1, "")
	if length(V2) > 0
		olExpPlotCurves!(p, V2, stroke2, showV2, ""; ls=:dash)
	end
	return p
end

# # Main plot
# plot(
# 	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")..., [120,160,190], readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [150,180,200]; title="Wing A1"), 
# 	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")..., [120,140,160], readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")..., [120,150,190]; title="Wing 4B1"),
# 	size=(800,400))

# plot of comparing different SDAB
p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.6), legend=:topleft, title="Different SDAB")
olExpPlotCurves!(p, readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [180], "Avik ")
olExpPlotCurves!(p, [180.0], [becky3L], [], "Becky 3L ")
olExpPlotCurves!(p, [180.0], [becky3R], [], "Becky 3R ")
olExpPlotCurves!(p, [180.0], [patrick], [], "Patrick ")
plot(p, size=(400,400))

gui()
	