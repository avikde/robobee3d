using Plots, DelimitedFiles

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlot(V, stroke; showV=[], title="")
	Nv = length(V)
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.55), legend=:topleft, title=title)
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=Int(V[i]), markershape=:auto, lw=2)
		end
	end
	return p
end

plot(
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 a1.csv")...; showV=[150,180,200], title="SDAB A1"),
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - halfbee1 4b1.csv")...; showV=[120,150,190], title="SDAB 4B1"),
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")...; showV=[120,160,190], title="Mod1 A1"), 
	olExpPlot(readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")...; showV=[120,140,160], title="Mod1 4B1"),
	size=(800,600))
gui()
	