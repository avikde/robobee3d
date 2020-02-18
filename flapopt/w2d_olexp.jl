using Plots, DelimitedFiles

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlot(V1, stroke1, showV1, V2, stroke2, showV2; title="")
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.6), legend=:topleft, title=title)

	function plotCurves(V, stroke, showV; kwargs...)
		Nv = length(V)
		for i=1:Nv
			if length(showV) == 0 || Int(V[i]) in showV
				plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=Int(V[i]), markershape=:auto, lw=2.5; kwargs...)
			end
		end
	end
	plotCurves(V1, stroke1, showV1)
	if length(V2) > 0
		plotCurves(V2, stroke2, showV2; ls=:dash)
	end
	return p
end

plot(
	# olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")...; showV=[150,180,200], title="SDAB A1"),
	# olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")...; showV=[120,150,190], title="SDAB 4B1"),
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")..., [120,160,190], readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")..., [150,180,200]; title="Wing A1"), 
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")..., [120,140,160], readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")..., [120,150,190]; title="Wing 4B1"),
	size=(800,400))
gui()
	