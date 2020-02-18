using Plots, DelimitedFiles

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlot(V, stroke; showV=[])
	Nv = length(V)
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.7), legend=:topleft)
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=Int(V[i]), markershape=:auto, lw=2)
		end
	end
	return p
end

mod1a1_V, mod1a1 = readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")
mod4bh1_V, mod4bh1 = readOLExpCSV("data/normstroke/Param opt manuf 2 - mod4 b h1.csv")

plot(
	olExpPlot(mod1a1_V, mod1a1; showV=[120,160,190]), 
	olExpPlot(mod4bh1_V, mod4bh1; showV=[120,150,160,170,180,190]))
gui()
	